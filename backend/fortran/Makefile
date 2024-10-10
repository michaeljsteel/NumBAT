win:
	meson setup --reconfigure --buildtype=release mbuild_win --native-file nb-win-native-file.ini
	meson compile -C mbuild_win
	copy mbuild_win\nb_fortran.dll nb_fortran.pyd



gccdebug:
	CC=gcc meson setup --buildtype=debugoptimized mbuild_gccdebug
	meson configure mbuild_gccdebug
	meson compile -C mbuild_gccdebug
	cp build_gccdebug/nb_fortran.so .

inteldebug:
	CC=icx meson setup --buildtype=debugoptimized mbuild_inteldebug
	meson compile -C mbuild_inteldebug

gcc:
	- rm -rf build_gcc
	meson setup --reconfigure --buildtype=release build_gcc 
	meson compile -v -C build_gcc
	cp build_gcc/nb_fortran.so .

intel:
	- rm -rf build_intel
	meson setup --reconfigure --buildtype=release build_intel \
		--native-file nb-linuxintel-native-file.ini 

	meson compile -C build_intel
	cp build_intel/nb_fortran.so .

mac: 
	- rm -rf build_mac
	meson setup --reconfigure --buildtype=release build_mac \
		--native-file nb-darwin-native-file.ini 

	meson compile -C build_mac
	cp build_mac/nb_fortran.dylib nb_fortran.so 

clean:
	-rm -rf build_*
	-rm nb_fortran.*
