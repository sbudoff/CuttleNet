import subprocess
import os

os.chdir('/home/sam/Poleg-Polsky/ICML/CellTICS')
 
for seed in range(18, 108, 18):
	# Format the dataset name and file paths based on the current seed
	dataset_name = f'Retina{seed}_Ablated'
	data_path = f'Retina/Retina_shuffle_ablated_{seed}/'
	reference_data_path = data_path + 'retina_rdata.csv'
	query_data_path = data_path + 'retina_qdata.csv'
	reference_label_path = data_path + 'retina_rlabel.csv'
 
	# Construct the command
	command = [
	    'python', '-u', 'code/main.py',
	    '-dataset_name', dataset_name,
	    '-reference_data_path', reference_data_path,
	    '-query_data_path', query_data_path,
	    '-reference_label_path', reference_label_path,
	    '-ensembl_pathway_relation', 'reactome/Ensembl2Reactome_All_Levels.txt',
	    '-pathway_names', 'reactome/ReactomePathways.txt',
	    '-pathway_relation', 'reactome/ReactomePathwaysRelation.txt',
	    '-n_hidden_layer', '5'
	]

	# Execute the command
	subprocess.run(command)
