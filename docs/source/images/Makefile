
.PHONY: images

all: meshfigs
	update_imgs.sh

meshfigs:
	python make_meshfigs.py
	- mv *-annotated.png *-refractive_index.png meshes
	- rm *.png

