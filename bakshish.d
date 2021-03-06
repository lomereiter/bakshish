/*
  Input BAM file is expected to contain reads with MD tags,
  mapped to the unique reference sequence.
  No 'N's or other ambiguous bases in reads are allowed.
  
  Memory usage: O(GENOMELENGTH)
*/

import bio.bam.reader, bio.bam.read, bio.bam.md.core,
    bio.core.base, bio.core.sequence, bio.core.fasta;
import std.exception, std.range, std.algorithm, std.stdio, std.typecons;

// Consider two mapped reads (forward strand):
// ...ACGTTT*TGTT...
// ...ACG*TTTTGTT...
// Common sense suggests that in both cases it should be
// ...ACGTTTT*GTT... instead, otherwise it's meaningless
// to collect any statistics about motifs preceding errors.
// The same goes for insertions inside homopolymers.
//
// Thus, we need to correct the alignment in such cases:
// 1) for each indel adjacent to two matching chunks, get the sequence;
// 2) if it's not a single homopolymer, skip the indel;
// 3) count how many nucleotides in the next aligned read chunk are
//    equal to the nucleotide in the inserted/deleted homopolymer;
// 4) increase length of the operation to the left by that amount;
// 5) decrease length of the operation to the right by that amount.
//
// Input: CIGAR to be fixed, query sequence, and list of deletions
// (each deletion is of type bio.core.sequence.NucleotideSequence)
//
// Note: deletions are retrieved from MD tags, it's assumed that
// the tags are correct and agree with CIGAR.
void fixCigar(Seq, Deletions)(ref CigarOperation[] cigar,
                              Seq nucleotides, Deletions deletions)
{
    if (cigar.length > 0) {
        if (cigar[0].is_query_consuming) nucleotides.popFrontN(cigar[0].length);
        if (cigar[0].type == 'D') deletions.popFront();
    }

    static bool singleNucIndel(S)(S seq) {
        auto nuc1 = seq[0];
        foreach (nuc2; seq[1 .. seq.length])
            if (nuc2 != nuc1)
                return false;
        return true;
    }

    for (size_t i = 1; i < cigar.length - 1; i++) {
        if ((cigar[i-1].type == '=' && cigar[i+1].type == '=') ||
            (cigar[i-1].type == 'M' && cigar[i+1].type == 'M')) {
            Base5 nuc;
            bool single_nuc_indel = false;
            if (cigar[i].type == 'I') {
                auto seq = nucleotides[0 .. cigar[i].length];
                nuc = cast(Base5)seq[0];
                single_nuc_indel = singleNucIndel(seq);
            } else if (cigar[i].type == 'D') {
                auto seq = deletions.front;
                nuc = cast(Base5)seq[0];
                single_nuc_indel = singleNucIndel(seq);
            }

            if (single_nuc_indel) {
                auto start = cigar[i].type == 'I' ? cigar[i].length : 0;
                auto end = start + cigar[i + 1].length;
                int change;
                foreach (nuc2; nucleotides[start .. end])
                    if (nuc2 == nuc) change += 1; else break;
                cigar[i - 1] = CigarOperation(cigar[i - 1].length + change,
                                              cigar[i - 1].type);
                cigar[i + 1] = CigarOperation(cigar[i + 1].length - change,
                                              cigar[i + 1].type);
                nucleotides.popFrontN(change);
            }
        }

        if (cigar[i].is_query_consuming)
            nucleotides.popFrontN(cigar[i].length);
        if (cigar[i].type == 'D')
            deletions.popFront();
    }
}

unittest {
    import std.stdio;
    writeln("Testing fixCigar...");
    CigarOperation[] ops = [CigarOperation(5, '='),
                            CigarOperation(2, 'I'),
                            CigarOperation(3, '=')];
    auto seq = nucleotideSequence("ACGGTTTTGA");
    //                                  **      before
    //                                   **     should be 
    NucleotideSequence[] dels;
    fixCigar(ops, seq, dels);
    assert(ops[0].length == 6);
    assert(ops[2].length == 2);

    seq = nucleotideSequence("ACGT" "TTGA" "A" "AAA");
    //                deletion of TT here;  ^-insertion   
    //
    //            should be: "ACGTTT" "GAAAA" "A" ""
    ops = [CigarOperation(4, '='), CigarOperation(2, 'D'),
           CigarOperation(4, '='), CigarOperation(1, 'I'),
           CigarOperation(3, '=')];
    dels = [nucleotideSequence("TT")];
    fixCigar(ops, seq, dels);
    assert(ops[0].length == 6);
    assert(ops[2].length == 5);
    assert(ops[4].length == 0);
}

// We count how many reads matched at each position,
// how many reads have mismatches at each position,
// how many deletions start at each position,
// and how many insertions start after each position.
struct SiteStats {
    uint insertions;
    uint deletions;
    uint mismatches;
    uint matches;

    uint errors() @property const {
        return insertions + deletions + mismatches;
    }

    void opOpAssign(string op)(ref const SiteStats other) if (op == "+") {
        insertions += other.insertions;
        deletions += other.deletions;
        mismatches += other.mismatches;
        matches += other.matches;
    }
}

alias StrandStats = SiteStats[];

struct GenomeStats {
    StrandStats forward;
    StrandStats reverse;

    this(size_t genome_length) {
        forward.length = genome_length;
        reverse.length = genome_length;
    }
}

CigarOperation[] localCigar(size_t N, Cigar)(ref CigarOperation[N] cigar_buf,
                                             ref size_t cigar_len,
                                             Cigar operations)
{
    auto ops = operations.save;
    while (!operations.empty) {
        cigar_buf[cigar_len++] = operations.front;
        operations.popFront();
        if (cigar_len == cigar_buf.length)
            return ops.map!(x => cast()x).array();
    }
    return cigar_buf[0 .. cigar_len];
}

void addRead(ref GenomeStats stats, BamRead read, bool flip) {
    assert(!read.is_unmapped);

    CigarOperation[256] cigar_buf = void;
    size_t cigar_len;
    auto cigar = localCigar(cigar_buf, cigar_len, read.extended_cigar);
    
    auto md = read["MD"];
    enforce(md.is_string);
    string md_str = *cast(string*)(&md);    
    auto deletions = mdOperations(md_str).filterBidirectional!q{ a.is_deletion }
                                         .map!q{ a.deletion };

    SiteStats[] strand_stats;
    size_t position;

    if (!flip) {
        position = read.position;
        fixCigar(cigar, read.sequence, deletions);
        strand_stats = read.strand == '+' ? stats.forward[] : stats.reverse[];
    } else {
        position = stats.forward.length - 1 - (read.position + read.basesCovered());
        reverse(cigar);
        fixCigar(cigar, read.sequence.retro.map!complementBase, deletions.retro);
        strand_stats = read.strand == '-' ? stats.forward[] : stats.reverse[];
    }

    while (!cigar.front.is_reference_consuming)
        cigar.popFront();

    foreach (op; cigar) {
        if (position >= strand_stats.length)
            break;
        
        switch (op.type) {
        case '=':
            foreach (i; 0 .. op.length) {
                if (position + i >= strand_stats.length)
                    break;
                strand_stats[position + i].matches += 1;
            }
            break;
        case 'X':
            foreach (i; 0 .. op.length) {
                if (position + i >= strand_stats.length)
                    break;
                strand_stats[position + i].mismatches += 1;
            }
            break;
        case 'D': strand_stats[position].deletions += 1; break;
        case 'I': strand_stats[position-1].insertions += 1; break;
        default: break;
        }

        if (op.is_reference_consuming) position += op.length;
    }
}

import dstats.tests;
import std.array, std.typecons;

void shift(ref uint kmer, char nuc, ubyte n) {
    kmer <<= 2;
    kmer &= (1 << (2*n)) - 1;
    kmer += Base5(nuc).internal_code;
}

auto advance(ref uint kmer, string genome, ref size_t pos,
             ref GenomeStats genome_stats, ubyte n) {
    SiteStats fwd, rev;
    while (true) {
        kmer.shift(genome[pos], n);
        fwd += genome_stats.forward[pos];
        rev += genome_stats.reverse[pos];
        ++pos;
        if (pos >= genome.length) break;
        if (genome[pos] != genome[pos-1]) break;
    }
    return Tuple!(SiteStats, "forward", SiteStats, "reverse")(fwd, rev);
}

void fillKmerStats(ref SiteStats[] kmer_fwd_stats,
                   ref SiteStats[] kmer_rev_stats,
                   ref GenomeStats stats, string genome, ubyte n)
{
    uint kmer = 0;
    size_t pos = 0;

    for (; pos < n - 1; ++pos) shift(kmer, genome[pos], n);
    
    while (true) {
        auto hstats = advance(kmer, genome, pos, stats, n);
        kmer_fwd_stats[kmer] += hstats.forward;
        kmer_rev_stats[kmer] += hstats.reverse;
        if (pos >= genome.length) break;
    }
}

static enum help_sections =
    ["main" :
     "Usage: bakshish [subcommand] [subcommand-specific options]\n"
     "       where subcommand is one of 'collect', 'merge', 'aggregate'\n"
     "\n"
     "       To get help about each subcommand, call it without options.\n"
     "\n"
     "       'collect' subcommand collects statistics about errors\n"
     "                 at the end of each k-mer\n"
     "\n"
     "       'merge' subcommand merges multiple lists produced by 'collect',\n"
     "               allowing to get more statistical power by considering\n"
     "               several BAM files\n"
     "\n"
     "       'aggregate' subcommand aims at searching motifs of length k\n"
     "                   with m 'N's inside, aggregating information from all\n"
     "                   k-mers matching each such pattern",

     "collect" :
     "Usage: bakshish collect <input.bam> <genome.fasta> <K>\n"
     "\n"
     "       Input BAM file must contain reads mapped to the same chromosome,\n"
     "       which is the first entry in the provided FASTA file; these reads\n"
     "       must also contain MD tags (TMAP writes them by default).\n"
     "\n"
     "       For each K-mer seen in the genome (and its reverse-complement)\n"
     "       the tool assesses strand bias in the errors at its end.\n"
     "\n"
     "       The results are output to <input.bam.kmer_stats.K> file\n"
     "       in table form."
     ];

void printHelp(string section="main") {
    stderr.writeln(help_sections[section]);
}

int main(string[] args) {
    import std.parallelism;
    defaultPoolThreads = 15;

    if (args.length == 0) {
        printHelp();
        return 0;
    }

    switch (args[1]) {
    case "collect":
        args = args[1 .. $];
        if (args.length < 4) {
            printHelp("collect");
            return -1;
        }
        auto n = args[3].to!ubyte;
        collectKmerStats(args[1], args[2], n, args[1] ~ ".kmer_stats." ~ n.to!string);
        break;
    case "merge":
    case "aggregate":
        stderr.writeln("NOT YET IMPLEMENTED");
        break;
    default:
        printHelp();
        break;
    }
    return 0;
}

double computePValue(uint kmer,
                     SiteStats[] kmer_fwd_stats, SiteStats[] kmer_rev_stats)
{
    uint[2][2] table = void;
    table[0][0] = kmer_fwd_stats[kmer].errors;
    table[0][1] = kmer_fwd_stats[kmer].matches;
    table[1][0] = kmer_rev_stats[kmer].errors;
    table[1][1] = kmer_rev_stats[kmer].matches;

    if (table[0][].chain(table[1][]).all!(x => x == 0))
        return 1;

    bool use_chisq = table[0][].chain(table[1][]).any!(x => x > 5000);

    double p;
    if (use_chisq) {
        uint[][2] t = void; t[0] = table[0][]; t[1] = table[1][];
        return chiSquareContingency(t[]).p;
    } else {
        return fisherExact(table).p;
    }
}

void collectKmerStats(string bam_filename, string fasta_filename, ubyte n,
                      string output_filename)
{
    auto bam = new BamReader(bam_filename);
    bam.assumeSequentialProcessing();
    auto genome = fastaRecords(fasta_filename).front.sequence;
    auto rgenome = genome.retro.map!(nuc => Base5(nuc).complement).array();

    auto output_file = File(output_filename, "w+");

    enforce(n >= 3, "k-mer length must be at least 3");
    enforce(n <= 12, "k-mer length must be at most 12");
    auto stats = GenomeStats(genome.length);
    auto rstats = GenomeStats(rgenome.length);

    foreach (read; bam.reads().filter!(r => !r.is_unmapped &&
                                            !r.is_secondary_alignment &&
                                       !r.is_duplicate))
    {
        stats.addRead(read, false);
        rstats.addRead(read, true);
    }

    auto kmer_fwd_stats = new SiteStats[](4 ^^ n);
    auto kmer_rev_stats = new SiteStats[](4 ^^ n);

    // collect per-kmer stats
    fillKmerStats(kmer_fwd_stats, kmer_rev_stats, stats, genome, n);
    fillKmerStats(kmer_fwd_stats, kmer_rev_stats, rstats, genome, n);
        
    output_file.writeln("kmer\tfm\tfmm\tfi\tfd\trm\trmm\tri\trd\tp");
    auto kmer_str = new char[n];
    foreach (kmer; 0 .. 4 ^^ n) {
        auto p = computePValue(kmer, kmer_fwd_stats, kmer_rev_stats);
        foreach (j; 0 .. n)
            kmer_str[n - j - 1] = Base5.fromInternalCode((kmer >> (2 * j)) & 3);

        output_file.writeln(kmer_str, "\t",
                            kmer_fwd_stats[kmer].matches, "\t",
                            kmer_fwd_stats[kmer].mismatches, "\t",
                            kmer_fwd_stats[kmer].insertions, "\t",
                            kmer_fwd_stats[kmer].deletions, "\t",
                            kmer_rev_stats[kmer].matches, "\t",
                            kmer_rev_stats[kmer].mismatches, "\t",
                            kmer_rev_stats[kmer].insertions, "\t",
                            kmer_rev_stats[kmer].deletions, "\t",
                            p);
    }
}
