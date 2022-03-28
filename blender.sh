source activate parl
tmux new -s parl
tmux attach -t parl
cd ~/ParlAI; python setup.py develop
python projects/personachat/scripts/kvmemnn_interactive.py
python parlai/scripts/safe_interactive.py -mf zoo:blender/blender_9B/model -t blended_skill_talk

