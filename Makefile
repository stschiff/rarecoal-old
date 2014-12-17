debug: build/debug/rarecoal

release: build/release/rarecoal

all: debug release

unittest: build/debug/unittest
	build/debug/unittest

build/debug/rarecoal : *.d
	dmd -debug $^ -odbuild/debug -of$@

build/release/rarecoal : *.d
	dmd -O -release $^ -odbuild/release -of$@

build/debug/unittest : *.d
	dmd -unittest $^ -odbuild/debug -ofbuild/debug/unittest

clean:
	rm -rf build/debug/* build/release/*

test: test.d model.d coal_state.d
	dmd test.d model.d coal_state.d

.PHONY: clean unittest