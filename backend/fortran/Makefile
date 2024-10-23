#UNAME_S := $(shell uname -s)


def:
	$(error 'Specify a target that is one of win, gcc, gccdebug, intel, inteldebug')



win:
#ifneq ($(UNAME_S), Windows)  # TODO: make this a make function?
#	$(error 'We do not seem to be building on Linux. Check you are choosing the right target')
#endif
	meson setup --reconfigure --buildtype=release build_win --native-file nb-win-native-file.ini
	meson compile -C build_win
	copy build_win\nb_fortran.dll nb_fortran.pyd



gccdebug:
	- rm -rf build_gccdebug
	CC=gcc meson setup --buildtype=debugoptimized build_gccdebug
	#meson configure build_gccdebug
	meson compile -C build_gccdebug
	cp build_gccdebug/nb_fortran.so .

inteldebug:
	CC=icx meson setup --buildtype=debugoptimized build_inteldebug
	meson compile -C build_inteldebug

gcc:
#ifneq ($(UNAME_S), Linux)  # TODO: make this a make function?
#	$(error 'We do not seem to be building on Linux. Check you are choosing the right target')
#endif
	- rm -rf build_gcc
	meson setup --reconfigure --buildtype=release build_gcc
	#meson compile -v -C build_gcc --jobs 1
	meson compile -v -C build_gcc
	cp build_gcc/nb_fortran.so .

intel:
	- rm -rf build_intel
	meson setup --reconfigure --buildtype=release build_intel \
		--native-file nb-linuxintel-native-file.ini

	meson compile -C build_intel
	cp build_intel/nb_fortran.so .

mac:
ifneq ($(UNAME_S), Darwin)  # TODO: make this a make function?
	$(error 'We do not seem to be building on Mac. Check you are choosing the right target')
endif
	- rm -rf build_mac
	meson setup --reconfigure --buildtype=release build_mac \
		--native-file nb-darwin-native-file.ini

	meson compile -C build_mac
	cp build_mac/nb_fortran.dylib nb_fortran.so

clean:
	-rm -rf build_*
	-rm nb_fortran.*
