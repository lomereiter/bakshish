debug:
	rdmd --force --build-only -IBioD bakshish.d

release:
	rdmd --force --build-only --compiler=ldmd2 -IBioD -O -release -inline bakshish.d
