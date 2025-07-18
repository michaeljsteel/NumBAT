NUMBAT_DIR = $HOME/numbat
NUMBAT_VER = latest

sudo apt-get update
sudo apt-get upgrade
sudo apt-get install gcc gfortran make gmsh python3-venv python3-dev meson pkg-
˓→config ninja-build

sudo add-apt-repository universe
sudo apt-get install libarpack2-dev libparpack2-dev libatlas-base-dev libblas-
˓→dev liblapack-dev libsuitesparse-dev

pip3 install numpy matplotlib scipy psutils

cd $HOME
python3 -m venv nbpy3
source ~/nbpy3/bin/activate

cd $NUMBAT_DIR

git clone https://github.com/michaeljsteel/NumBAT.git $NUMBAT_VER
cd $NUMBAT_VER

cd backend/fortran

make gcc

cd ..
python ./nb_install_tester.py

echo "NumBAT installation complete."