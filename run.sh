source venv/bin/activate
ProTECT --config mustard_config.yaml --workDir /scratch/drkthomp/workDir /scratch/drkthomp/jobStore --restart|& tee errors/$(date '+%Y-%m-%d-%H-%M-%S').txt

