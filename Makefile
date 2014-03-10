debug:
	rdmd --force --build-only -IBioD nghytser.d

release:
	rdmd --force --build-only --compiler=ldmd2 -IBioD -O -release -inline nghytser.d
