rarecoal : *.d
	dmd -O $^ -of$@

unittest : *.d
	rdmd -unittest main.d
