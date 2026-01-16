# This is the single point of definition of NumBAT version numbers

NUMBAT_VERSION_MAJ=2
NUMBAT_VERSION_MIN=1
NUMBAT_VERSION_MIC=5

######################################
# No changes below here

NUMBAT_VERSION_STR_MMM = f'{NUMBAT_VERSION_MAJ}.{NUMBAT_VERSION_MIN}.{NUMBAT_VERSION_MIC}'
NUMBAT_VERSION_STR_MM = f'{NUMBAT_VERSION_MAJ}.{NUMBAT_VERSION_MIN}'


if __name__ == '__main__':

    with open('src/nbversion_incl.h', 'wt') as fout:
        fout.write(f'character(len=100), parameter :: NUMBAT_VERSION_STR_MMM = "{NUMBAT_VERSION_STR_MMM}"\n')
        fout.write(f'character(len=100), parameter :: NUMBAT_VERSION_STR_MM = "{NUMBAT_VERSION_STR_MM}"\n')
        fout.write(f'integer(8), parameter :: NUMBAT_VERSION_MAJ = {NUMBAT_VERSION_MAJ}\n')
        fout.write(f'integer(8), parameter :: NUMBAT_VERSION_MIN = {NUMBAT_VERSION_MIN}\n')
        fout.write(f'integer(8), parameter :: NUMBAT_VERSION_MIC = {NUMBAT_VERSION_MIC}\n')




    print(NUMBAT_VERSION_STR_MMM) # for meson to pick up
