debug: build/debug/rarecoal

release: build/release/rarecoal

unittest:
	build/debug/unittest

build/debug/rarecoal : *.d
	dmd -debug $^ -odbuild/debug -of$@

build/release/rarecoal : *.d
	dmd -O -release $^ -odbuild/release -of$@

build/debug/unittest : *.d
	dmd -unittest $^ -odbuild/debug -ofbuild/debug/unittest

clean:
	rm -rf build/debug/* build/release/*

.PHONY: clean unittest