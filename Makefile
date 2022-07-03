.PHONY: env
env: 
	mamba env create -f environment.yml -p ~/envs/dev_env
	conda activate /home/jovyan/envs/home/jovyan/hw07-hw07-group12/dev_env;python -m ipykernel install --user --name=dev_env

.PHONY : clean
clean :
	rm -f figures/*.png

.PHONY : all
all :
	jupyter execute RNA_velocity_concept.ipynb
	jupyter execute main.ipynb
	jupyter execute dentate_vignette.ipynb
	jupyter execute dentategyrus_JQ.ipynb