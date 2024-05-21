#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pyreadr, pickle, csv, re, os, tqdm
import numpy as np
import pandas as pd
import tifffile

import torch
import torch.nn as nn

from sklearn.preprocessing import LabelEncoder
from sklearn.utils import shuffle

import matplotlib.pyplot as plt


# In[2]:


# Data Dependencies
path = '/home/sam/scRNAseq/Xenium/Network_genes_NoiseInjection.RData'
model_path = "/home/sam/scRNAseq/Xenium/AlonNN/NoiseInj/model_state_epochs_150_earlyStop_50_l1_0.0001_depth_5_withSkips_seed_18.pt"


# # Model Loading And Configuration

# In[3]:


# Load scRNAseq Dataset 
rdata = pyreadr.read_r(path)
# Load data
df = rdata['Retina_expMatrix_candidateGenes']
df['Cluster'] = df['Cluster'].apply(lambda x: x if len(x.split('_')[0]) == 2 else '0' + x) # Standardize cluster names


# Load the list of indices for each network to use
class_net_genes = rdata['Class_indices'].to_numpy().ravel()
rgc_net_genes = rdata['RGC_indices'].to_numpy().ravel()
ac_net_genes = rdata['AC_indices'].to_numpy().ravel()
bc_net_genes = rdata['BC_indices'].to_numpy().ravel()
hc_net_genes = rdata['HC_indices'].to_numpy().ravel()
nn_net_genes = rdata['NonNeural_indices'].to_numpy().ravel()

# Encode and format scRNAseq dat
df['cluster'] = df['Cluster']

def encode_class(arr):
    '''This function will encode subtypes' cell classesbased on expert rules and is not intended for decoding'''
    custom_array = []

    for value in arr:
        if re.match(r'^\d{2}_', value):
            custom_array.append(0)
        elif value.startswith('AC_'):
            custom_array.append(1)
        elif value.endswith('Photoreceptors'):
            custom_array.append(2)
        elif value == '0MG (Mueller Glia)':
            custom_array.append(3)
        elif value.startswith('0BC'):
            custom_array.append(4)
        elif value.startswith('0RBC'):
            # Note this duplication is for simplicity of handling the 2 BC naming conventions
            custom_array.append(4)
        elif value == '0Horizontal Cell':
            custom_array.append(5)
        elif value == '0Pericyte':
            custom_array.append(6)
        elif value == '0Endothelial':
            custom_array.append(6)
        elif value == '0Microglia':
            custom_array.append(6)   
        else:
            custom_array.append(7)
    return custom_array

# Function to generate a sort key for each subclass name
def sort_key(name):
    if name.startswith('0BC'):
        return 4
    elif name.startswith('0RBC'):
        return 4
    elif re.match(r'^\d{2}_', name):
        return 0
    elif name.startswith('AC_'):
        return 1
    elif name.endswith('Photoreceptors'):
        return 2
    elif name == '0MG (Mueller Glia)':
        return 3
    elif name == '0Horizontal Cell':
        return 5
    elif name == '0Pericyte':
        return 6
    elif name == '0Endothelial':
        return 6
    elif name == '0Microglia':
        return 6
    else:
        return 7

# Apply the sort_key function to each subclass name and sort the DataFrame
df['sort_key'] = df['Cluster'].apply(sort_key)
df.sort_values(by='sort_key', inplace=True)
df.drop(columns='sort_key', inplace=True)  # Optionally remove the sort key

class_arr = encode_class(df['Cluster'])

# Encode the categoric response 
le = LabelEncoder()
df['Cluster'] = le.fit_transform(df['Cluster'])

cluster_col = df.pop('Cluster')
dataset_col = df.pop('Dataset')
df.insert(len(df.columns), 'Cluster', cluster_col)
df.insert(len(df.columns), 'Class', class_arr)

def create_mapping(df):
    # Extract the unique pairs of encoded cluster values and their corresponding class encodings
    unique_pairs = df[['Cluster', 'Class']].drop_duplicates()
    
    # Create a dictionary mapping from Cluster to Class
    mapping = dict(zip(unique_pairs['Cluster'], unique_pairs['Class']))
    
    return mapping

# Usage
mapping = create_mapping(df)

print(mapping)


# Utilize the mapping to count the number of subclasses for each class
subclass_counts = {i: 0 for i in range(8)}  # Initialize counts for 6 classes
for _, class_id in mapping.items():
    subclass_counts[class_id] += 1

# Construct the class_info dictionary
class_info = {
    'Genes': class_net_genes,  # Genes used for class classification
    'num_classes': 8,  # Total number of classes
    0: {  # Information for class 0 (RGCs)
        'Genes': rgc_net_genes,
        'num_subclasses': subclass_counts[0]
    },
    1: {  # Information for class 1 (ACs)
        'Genes': ac_net_genes,
        'num_subclasses': subclass_counts[1]
    },
    2: {  # Information for class 2
        'Genes': bc_net_genes,
        'num_subclasses': subclass_counts[2]
    },
    3: {  # Information for class 3
        'Genes': bc_net_genes,
        'num_subclasses': subclass_counts[3]
    },
    4: {  # Information for class 4BC
        'Genes': bc_net_genes,
        'num_subclasses': subclass_counts[4]
    },
    5: {  # Information for class 5 HCs
        'Genes': hc_net_genes,
        'num_subclasses': subclass_counts[5]  
    },
    6: {  # Information for class 6 NonNeural
        'Genes': nn_net_genes,
        'num_subclasses': subclass_counts[6]  
    },
    7: {  # Information for class 5 (Catch-all class)
        'Genes': class_net_genes,
        'num_subclasses': 1  # Only 1 subclass as it's a catch-all class
    }
}

# Add 'num_hidden' with default zero to each tentacle
for c in range(class_info['num_classes']):
    class_info[c]['num_hidden'] = 0
    class_info[c]['skip'] = False
    

l_arm = 5
skip = True
class_info[0]['num_hidden'] = l_arm # RGC Arm
class_info[1]['num_hidden'] = l_arm # AC Arm
class_info[0]['skip'] = skip # RGC Arm
class_info[1]['skip'] = skip # AC Arm


# In[4]:


class TentacleNet(nn.Module):
    def __init__(self, input_size, num_subclasses, num_hidden, skip = False):
        super(TentacleNet, self).__init__()
        self.fc1 = nn.Linear(input_size, 2*num_subclasses)
        self.num_hidden = num_hidden
        self.skip = skip+0
        if self.num_hidden > 0:
            self.hidden = nn.ModuleList([nn.Linear(2*num_subclasses, 2*num_subclasses) for _ in range(num_hidden)])
        self.fc2 = nn.Linear(2*num_subclasses, num_subclasses)

    def forward(self, x):
        x = nn.functional.relu(self.fc1(x))
        if self.num_hidden > 0:
            x_skip = x*self.skip  # Save output of fc1 for skip connection
            for hidden_layer in self.hidden:
                x = nn.functional.relu(hidden_layer(x))
            x = x + x_skip  # Add skip connection before final activation
        x = self.fc2(x)
        return nn.functional.log_softmax(x, dim=1)


class CuttleNet(nn.Module):
    def __init__(self, class_info, mapping):
        super(CuttleNet, self).__init__()
        
        self.device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        self.n = len(mapping) # Number of subclasses
        self.class_info = class_info

        # Class Classifier
        self.class_fc1 = nn.Linear(len(class_info['Genes']), 2*class_info['num_classes'])
        self.class_fc2 = nn.Linear(2*class_info['num_classes'], class_info['num_classes'])

        # Subclass Classifiers
        self.subclass_nets = nn.ModuleDict({
            str(class_id): TentacleNet(input_size=len(subclass_info['Genes']) + class_info['num_classes'], 
                                       num_subclasses=subclass_info['num_subclasses'],
                                      num_hidden = subclass_info['num_hidden'],
                                      skip = subclass_info['skip'])
            for class_id, subclass_info in class_info.items()
            if isinstance(class_id, int)
        })
        
        # Calculate the number of subclasses for each class
        self.num_subclasses_per_class = self.calculate_subclasses_per_class(mapping)
        
    def get_subclass_range_for_class(self, class_id):
        start_index = sum(self.num_subclasses_per_class[cid] for cid in range(class_id))
        end_index = start_index + self.num_subclasses_per_class[class_id]
        return slice(start_index, end_index)
    
    def calculate_subclasses_per_class(self, mapping):
        """
        Calculate the number of subclasses for each class using the mapping.
        """
        num_subclasses_per_class = {class_id: 0 for class_id in range(self.class_info['num_classes'])}
        for subclass_id in mapping.keys():
            class_id = mapping[subclass_id]
            num_subclasses_per_class[class_id] += 1
        return num_subclasses_per_class

    def forward(self, x):
        # Class classification
        class_genes = x[:, self.class_info['Genes']]
        class_x = nn.functional.relu(self.class_fc1(class_genes))
        class_output = nn.functional.log_softmax(self.class_fc2(class_x), dim=1)

        # Initialize an output tensor for all subclasses
        all_subclass_output = torch.zeros(x.size(0), 130, device=self.device)  # Assuming 130 total subclasses

        # Populate the output tensor
        for class_id, subclass_info in self.class_info.items():
            if isinstance(class_id, int):
                subclass_genes = x[:, subclass_info['Genes']]
                subclass_input = torch.cat((subclass_genes, class_output), dim=1)

                # Convert class_id to string
                class_id_str = str(class_id)
                subclass_output = self.subclass_nets[class_id_str](subclass_input)

                # Get the range for this class's subclasses
                subclass_range = self.get_subclass_range_for_class(class_id)

                # Multiply subclass predictions by the class prediction probability
                all_subclass_output[:, subclass_range] = subclass_output * class_output[:, class_id].unsqueeze(1)

        return all_subclass_output

    
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

model = CuttleNet(class_info=class_info, mapping=mapping)

# Load the model state
model.load_state_dict(torch.load(model_path, map_location=device))

# Move the model to the appropriate device and set it to evaluation mode
model.to(device)
model.eval()


# # Inference

# In[5]:


# Step 1: Select the slice
df_slice = df.iloc[:, -3:-1]

# Step 2: Drop the index
df_slice_reset = df_slice.reset_index(drop=True)

# Step 3: Remove duplicate rows
df_unique = df_slice_reset.drop_duplicates()

# Step 4: Sort by the last column (you can reference it by its position since it's a slice)
df_sorted = df_unique.sort_values(by=df_unique.columns[-1])
clust_ids = list(df_sorted[['cluster']].values.ravel())
clust_ids.append('other')


# Extract Correct gene order
correct_order = np.array(df.iloc[:,0:-3].columns)
# Correct any mislabeled genes
gene_rename_map = {
'Gm11744': 'Prcd',
'Fam19a3': 'Tafa3',
#     'A730046J19Rik': 'Sertm2',
'Fam19a1': 'Tafa1',
'Cyr61': 'Ccn1'
}
correct_order = np.array([gene_rename_map.get(gene, gene) for gene in correct_order])


def create_count_matrix(df, correct_order):
    # Pivot table to count occurrences of feature_name for each cell_id
    count_matrix = pd.pivot_table(df, index='cell_id', columns='feature_name', aggfunc='size', fill_value=0)
    count_matrix = count_matrix.loc[:,correct_order] # extract genes only
    
    return count_matrix

def CuttleNet_Inference(data_root, clust_ids, correct_order,
                        chunk_size = 10000, visualize=False):


	print(f'Loading experiment data from {data_root}')
	# load transcript data
	data_path = os.path.join(data_root,"transcripts.csv.gz")
	# cell_shape_path = os.path.join(data_root,"cell_boundaries.csv.gz")
	# nuc_shape_path = os.path.join(data_root,"nucleus_boundaries.csv.gz")
	cent_shape_path = os.path.join(data_root,"cells.csv.gz")
	# Reading the compressed CSV file as chuncks
	# Create iterators for each dataset
	cs = 10**7
	# cell_shape_iter = pd.read_csv(cell_shape_path, compression='gzip', chunksize=cs)
	# nuc_shape_iter = pd.read_csv(nuc_shape_path, compression='gzip', chunksize=cs)
	# cent_shape_iter = pd.read_csv(cent_shape_path, compression='gzip', chunksize=cs)
	# Create an empty list to hold chunks
	data_chunks = []
	n_cells = 0
	for xen_data in pd.read_csv(data_path, compression='gzip', chunksize=cs):

		data = create_count_matrix(xen_data, correct_order)

		# Remove the row with the index 'UNASSIGNED'
		data = data.drop('UNASSIGNED')
		# Normalize each column by its maximum value, multiply by 1000, round, and divide by 1000
		data = data.apply(lambda x: round(1000 * x / x.max()) / 1000)
		# store cell ids
		cell_ids = list(data.index)

		# Convert to torch tensor and store on GPU
		expMatrix = data.to_numpy()
		expMatrix = torch.tensor(expMatrix, dtype=torch.float32)

		print(f'Data loaded containing {len(cell_ids)} cells')
		n_cells += len(cell_ids)

		# Calculate the number of chunks needed
		n_chunks = int(np.ceil(expMatrix.size(0) / chunk_size))

		# Placeholder to collect the output
		results = []

		print('Performing Inference')
		# Process each chunk
		for i in tqdm.tqdm(range(n_chunks)):
			# Calculate the start and end indices of the current chunk
			start_idx = i * chunk_size
			end_idx = min((i + 1) * chunk_size, expMatrix.size(0))

			# Extract the chunk
			chunk = expMatrix[start_idx:end_idx]

			# Move the chunk to GPU
			chunk = chunk.to('cuda')

			# Perform inference
			with torch.no_grad():  # Ensure gradients are not computed to save memory
			    chunk_output = model(chunk)

			# Move the results back to CPU and store them
			chunk_output = chunk_output.cpu()
			results.append(chunk_output)

		# Concatenate the results into a single tensor
		final_results = torch.cat(results, dim=0)

		final_df = pd.DataFrame(final_results.numpy(), columns= clust_ids)
		# Add Prediction column
		final_df['Prediction'] = final_df.idxmax(axis=1)
		# Add cell_ids column
		final_df['cell_id'] = cell_ids

		# Load and store cell shape information    
		# print('Loading shape information')

		# Load corresponding chunks for other datasets
		# cell_shape_data = next(cell_shape_iter)
		# nuc_shape_data = next(nuc_shape_iter)
		# cent_shape_data = next(cent_shape_iter)

		# Add suffix and rename columns as required
		# cell_shape_data = cell_shape_data.add_suffix('_cell').rename(index=str, columns={'cell_id_cell':'cell_id'})
		# nuc_shape_data = nuc_shape_data.add_suffix('_nucleus').rename(index=str, columns={'cell_id_nucleus':'cell_id'})
		

		# Merge data for shape information
		# bounds_shape_data = cell_shape_data.merge(nuc_shape_data, on='cell_id')
		# bounds_shape_data = bounds_shape_data.merge(cent_shape_data, on='cell_id')

		
		data_chunks.append(final_df)


	# Concatenate chunks into a single DataFrame
	print(f'Full dataset contains {n_cells} cells')
	final_df = pd.concat(data_chunks, axis=0)
	
	print("Loading Centroids")
	cent_shape = pd.read_csv(cent_shape_path, compression='gzip')
	cent_shape = cent_shape[['cell_id', 'x_centroid', 'y_centroid']]
	# Merge with the main dataframe chunk
	final_df = xen_data.merge(cent_shape, on='cell_id')


	print('Inference and dataframe merging complete.')

	return final_df


# In[ ]:


experiments = {0 : {'slide' : '0018429',
                   'path' : '/media/sam/New Volume/Xenium_Data/output-XETG00230__0018429__Region_1__20240105__233208'},
               1 : {'slide' : '0018432',
                   'path' : '/media/sam/New Volume/Xenium_Data/output-XETG00230__0018432__Region_2__20240105__233208'},
               2 : {'slide' : '0018336',
                   'path' : '/media/sam/New Volume/Xenium_Data/BudoffRun2_Slide 3_4/BudoffRun2_Slide 3_4/output-XETG00230__0018336__Region_1__20240124__002923'},
               3 : {'slide' : '0018521',
                   'path' : '/media/sam/New Volume/Xenium_Data/BudoffRun2_Slide 3_4/BudoffRun2_Slide 3_4/output-XETG00230__0018521__Region_1__20240124__002923'},
               4 : {'slide' : '0018624',
                   'path' : '/media/sam/New Volume/Xenium_Data/BudoffRun3_Slide 5_6/BudoffRun3_Slide 5_6/output-XETG00230__0018624__Region_1__20240127__000149'},
               5 : {'slide' : '0022826',
                   'path' : '/media/sam/New Volume/Xenium_Data/BudoffRun3_Slide 5_6/BudoffRun3_Slide 5_6/output-XETG00230__0022826__Region_1__20240127__000149'},
               6 : {'slide' : '0018300',
                   'path' : '/media/sam/New Volume/Xenium_Data/BudoffRun4_Slide 7_8/BudoffRun4_Slide 7_8/output-XETG00230__0018300__Region_1__20240206__235339'},
               7 : {'slide' : '0022825',
                   'path' : '/media/sam/New Volume/Xenium_Data/BudoffRun4_Slide 7_8/BudoffRun4_Slide 7_8/output-XETG00230__0022825__Region_1__20240206__235339'}}

save_path = '/home/sam/scRNAseq/Xenium/Full_Inference_All_Experiments.csv'

for i in range(len(experiments)):
    data_root = experiments[i]['path']
    inf_df = CuttleNet_Inference(data_root, clust_ids, correct_order)
    inf_df['slide'] = experiments[i]['slide']
    if i == 0:
        full_df = inf_df.copy()
    else:
        full_df = pd.concat((full_df, inf_df))
    
    print('Saving full dataframe')
    full_df.to_csv(save_path)

