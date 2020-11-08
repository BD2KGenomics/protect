source venv/bin/activate
sudo ProTECT --config ProTECT_config.yaml  --workDir /home/drkthomp/workDir /home/drkthomp/d/jobStore --restart |& tee errors/$(date '+%Y-%m-%d-%H-%M-%S').txt

