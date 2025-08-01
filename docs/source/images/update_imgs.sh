
mkdir -p tutorial
mkdir -p josab_tutorial
mkdir -p lit_examples

for tdir in  02 06a 06b 07 08 11a
    #14_0 14_5 14_10
do
mkdir -p tutorial/tut_${tdir}-fields
done

BASE="../../../../examples"

cd tutorial

for im in tut_01-mesh.png tut_01-refractive_index.png \
        tut_02-gain_spectra-logy.png tut_02-gain_spectra.png tut_02-gain_spectra_zoom.png  \
        tut_02-fields/AC_mode_03.png tut_02-fields/AC_mode_04.png tut_02-fields/EM_E_mode_00.png tut_02-fields/EM_H_mode_00.png \
        tut_03a-dispersion_symmetrised.png tut_03b-dispersion_multicore.png \
        tut_06a-gain_spectra.png tut_06a-gain_spectra-logy.png tut_06a-refractive_index.png tut_06b-refractive_index.png \
        tut_06a-fields/AC_mode_05.png tut_06a-fields/AC_mode_08.png tut_06a-fields/EM_E_mode_00.png \
        tut_08-acdisp_slotwidth.png \
        tut_11a-gain_spectra.png tut_11b-gain_spectra.png \
        tut_11a-fields/AC_mode_17.png tut_11a-fields/AC_mode_32.png tut_11a-fields/AC_mode_39.png tut_11a-fields/AC_mode_44.png tut_11a-fields/AC_mode_66.png tut_11a-fields/AC_mode_69.png tut_11a-fields/EM_E_mode_00.png  \
        tut_12-sil-emdisp_ng.png tut_12-sil-emdisp_Vb.png tut_12-sil-emdisp_Vneff.png tut_12-smf28-emdisp_ng.png tut_12-smf28-emdisp_Vb.png tut_12-smf28-emdisp_Vneff.png \
        tut_13-smf28-acdisp_qneff.png tut_13-smf28-acdisp_qnu.png

# tut_09-gain_spectra.png \
#        tut_14_0-fields tut_14_10-fields tut_14_5-fields tut_14-emdisp_ng.png tut_14-emdisp_Vneff.png \
#        tut_14_0-fields/EM_E_mode_00.png tut_14_0-fields/EM_E_mode_01.png tut_14_0-fields/EM_E_mode_02.png tut_14_0-fields/EM_E_mode_03.png tut_14_10-fields/EM_E_mode_00.png tut_14_10-fields/EM_E_mode_01.png tut_14_10-fields/EM_E_mode_02.png tut_14_10-fields/EM_E_mode_03.png tut_14_5-fields/EM_E_mode_00.png tut_14_5-fields/EM_E_mode_01.png tut_14_5-fields/EM_E_mode_02.png tut_14_5-fields/EM_E_mode_03.png
    do
        cp -R ${BASE}/tutorials/$im $im
    done

    mkdir -p tut_04-out tut_05-out tut_07-out tut_08-out

    for im in tut_04-out/tut_04-gain_spectra-scan.png \
        tut_05-out/tut_05-convergence-freq_AC.png tut_05-out/tut_05-convergence-freq_EM.png \
        tut_05-out/tut_05-convergence-gain_MB.png tut_05-out/tut_05-convergence-gain_PE.png \
        tut_05-out/tut_05-convergence-gain.png \
        tut_07-out/tut_07-gain_spectra.png tut_07-out/tut_07-gain_spectra-logy.png \
        tut_07-out/tut_07-refractive_index.png tut_07-out/tut_07-refractive_index_xcut.png \
        tut_07-out/tut_07-refractive_index_ycut.png tut_07-out/tut_07-elastic_velocity_v0.png \
        tut_07-out/tut_07-elastic_velocity_xcut.png  tut_07-out/tut_07-elastic_velocity_linecut.png
    do
        cp -pR ${BASE}/tutorials/$im $im
    done

    cp -pR ${BASE}/tutorials/tut_07-out/tut_07-fields/EM_E_mode_00.png tut_07-out
    cp -pR ${BASE}/tutorials/tut_07-out/tut_07-fields/AC_mode_00.png  tut_07-out
    cp -pR ${BASE}/tutorials/tut_07-out/tut_07-fields/AC_mode_05.png  tut_07-out
    cp -pR ${BASE}/tutorials/tut_07-out/tut_07-fields/AC_mode_05_ycut.png  tut_07-out
    cp -pR ${BASE}/tutorials/tut_07-out/tut_07-fields/AC_mode_12.png  tut_07-out
    cp -pR ${BASE}/tutorials/tut_08-out/tut_08_150.0-fields/AC_mode_00.png tut_08-out
    cp -pR ${BASE}/tutorials/tut_08-out/tut_08_150.0-fields/AC_mode_00_uy.png tut_08-out
    cp -pR ${BASE}/tutorials/tut_08-out/tut_08_150.0-fields/AC_mode_00_uz.png tut_08-out
    cp -pR ${BASE}/tutorials/tut_08-out/tut_08_150.0-fields/AC_mode_01.png tut_08-out
    cp -pR ${BASE}/tutorials/tut_08-out/tut_08_150.0-fields/AC_mode_02.png tut_08-out

cd ..

#exit

for ndir in  01 02a 03 04 05 06
do
mkdir -p josab_tutorial/josab_${ndir}-fields
done


cd josab_tutorial

for im in josab_01-fields/EM_E_mode_01.png josab_01-fields/AC_mode_28.png josab_01-gain_spectra.png \
    josab_02a-fields/EM_E_mode_00.png josab_02a-fields/AC_mode_02.png josab_02a-fields/AC_mode_06.png josab_02a-gain_spectra.png  \
    josab_02b-disp-qnu.png josab_02b-disp-qnu.png josab_02b-disp-qnu.png \
    josab_03-gain_spectra.png josab_03-fields/EM_E_mode_01.png josab_03-fields/AC_mode_06.png josab_03-fields/AC_mode_12.png \
    josab_04-gain_spectra.png josab_04-fields/EM_E_mode_00.png josab_04-fields/AC_mode_06.png \
    josab_05-gain_spectra.png josab_05-fields/EM_E_mode_00.png josab_05-fields/EM_E_mode_01.png josab_05-fields/AC_mode_06.png \
    josab_06-gain_spectra.png josab_06-fields/AC_mode_02.png josab_06-fields/EM_E_mode_00.png josab_06-fields/EM_E_mode_02.png
    do
        cp -pR ${BASE}/josab_tutorial/$im $im
    done

cd ..

for ndir in 01 02 03-diam-1050 04b 05 06a 07 08
do
mkdir -p lit_examples/lit_${ndir}-fields
done
mkdir -p lit_examples/lit_03-out/lit_03-diam-1050-fields

cd lit_examples
#lit_03-gain_spectra-logy_w1000.png lit_03-gain_tot-diam_scan.png \
for im in lit_01-gain_spectra-logy.png lit_01-fields/EM_E_mode_00.png \
    lit_01-fields/AC_mode_04.png lit_01-fields/AC_mode_55.png \
    lit_01-fields/AC_mode_52.png \
    lit_02-gain_spectra-logy.png lit_02-fields/AC_mode_23.png lit_02-fields/AC_mode_296.png \
    lit_03-gain_tot-diam_scan.png \
    lit_03-dwide-gain_tot-diam_scan.png \
    lit_03-exactdisp_wide-acdisp_qnu_exact.png \
    lit_03-out/lit_03-gain_spectra-logy_w1050.png \
    lit_03-out/lit_03-diam-1050-fields/EM_E_mode_00.png \
    lit_03-out/lit_03-diam-1050-fields/AC_mode_02.png \
    lit_03-out/lit_03-diam-1050-fields/AC_mode_05.png \
    lit_03-out/lit_03-diam-1050-fields/AC_mode_21.png \
    lit_03-out/lit_03-diam-1050-fields/AC_mode_28.png \
    lit_04b-fields/EM_E_mode_00.png \
    lit_04b-fields/AC_mode_38.png lit_04a-gain_spectra.png \
    lit_05-fields/EM_E_mode_00.png lit_05-fields/AC_mode_06.png \
    lit_06a-fields/AC_mode_04.png lit_06a-fields/AC_mode_05.png \
    lit_06a-gain_spectra.png lit_06a-gain_spectra-5.png \
    lit_06a-gain_spectra-6.png lit_06a-gain_spectra-8.png \
    lit_06a-gain_spectra-11.png lit_06b-gain_spectra.png \
    lit_06b-gain_spectra-logy.png lit_07-fields/EM_E_mode_00.png \
    lit_07-fields/AC_mode_19.png lit_07-gain_spectra.png \
    lit_08-fields/EM_E_mode_00.png lit_08-fields/EM_E_mode_01.png \
    lit_08-fields/AC_mode_23.png lit_08-gain_spectra.png \
    lit_09-gain_spectra.png lit_10a-gain_spectra.png lit_10b-gain_spectra.png
do
    cp -pR ${BASE}/lit_examples/$im $im
done
