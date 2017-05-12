.PHONY = all release clean test

all: release

release:
	mkdir -p _build && \
	cd _build && \
	cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 .. && \
	make

test: release
	cd _build/test && \
	./test_ptw

clean:
	rm -rf _build
