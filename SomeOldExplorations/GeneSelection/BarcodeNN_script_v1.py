import pyreadr, pickle
import numpy as np
import pandas as pd

import torch
import torch.nn as nn
# import torch.optim as optim
from torch.utils.data import DataLoader, random_split#, TensorDataset

from sklearn.preprocessing import LabelEncoder
from sklearn.utils import shuffle

############################################# HYPERPARAMETERS#####################################################################################
##################################################################################################################################################
# Inits
# Set the device to use (GPU if available, otherwise CPU)
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

num_epochs = 100 # specify the number of epochs to train for
batch_size = 32 # specify the batch size for training
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

curriculum_dict = {}

path = '/home/sam/Classes/Stats/Bayes/assignments/keyGeneExpressionBySubtype.RData'
path = '/home/sam/scRNAseq/RNAscope/RNAscope/RGC_df354.Rdata'

############################################# FUNCTIONS ##########################################################################################
##################################################################################################################################################
class Net(nn.Module):
    def __init__(self, input_size, num_classes,
                 h1_size):
        super(Net, self).__init__()
        self.fc1 = nn.Linear(input_size, h1_size)
        self.fc2 = nn.Linear(h1_size, num_classes)

    def forward(self, x):
        x = nn.functional.relu(self.fc1(x))
        x = self.fc2(x)
        return nn.functional.log_softmax(x, dim=1)
    
#################################################################################################################################################################################

def head2head(train_y, train_X, desired, noise):
    '''This function will take a desired category and a noise category and create a paired down dataset
    The paired down data set will be an equal mix of desired and noise
    If there are not enough noise examples, a random set will be duplicated to make the lengths match
    If there are too many noise samples, a random set will be excluded'''

    # Isolate the desired instances
    TS_desired_y = train_y[train_y == desired]
    TS_desired_X = train_X[train_y == desired]
    TS_noise_y = train_y[train_y == noise]
    TS_noise_X = train_X[train_y == noise]
    # Check if the desired class is longer, if so balance it with replicates from te noise
    if len(TS_desired_y) > len(TS_noise_y):
        delta = len(TS_desired_y) - len(TS_noise_y)
        indice = np.random.choice( range(len(TS_noise_y)), delta)
        indice = torch.tensor(indice)
        extra_y = TS_noise_y[indice]
        extra_X = TS_noise_X[indice]
        TS_noise_y = torch.cat((TS_noise_y,extra_y))
        TS_noise_X = torch.cat((TS_noise_X,extra_X))
    # Check if the desired class is shorter, if so randomely discard examples from the other
    elif len(TS_desired_y) < len(TS_noise_y):
        indice = np.random.choice( range(len(TS_noise_y)), len(TS_desired_y))
        indice = torch.tensor(indice)
        TS_noise_y = TS_noise_y[indice]
        TS_noise_X = TS_noise_X[indice]
    # Combine the two balanced data sets
    TS_y = torch.cat((TS_noise_y,TS_desired_y))
    TS_X = torch.cat((TS_noise_X,TS_desired_X))
    # Shuffle the data
    indice = np.random.choice( range(len(TS_y)), len(TS_y))
    indice = torch.tensor(indice)
    TS_y = TS_y[indice]
    TS_X = TS_X[indice]
    return TS_y, TS_X

#################################################################################################################################################################################

def DesiredGeneCurriculum(train_y, train_X, desired=1):
    if desired!=None:
        # Create a training set using this function to create my desired class biased training set
        noise_inds = [i for i in np.unique(train_y) if i != desired]

    # Move the data to the GPU
    # Set the device to use (GPU if available, otherwise CPU)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    # Convert the input data to the appropriate data type for the GPU
    train_X = train_X.to(device)
    train_y = train_y.to(device)

    ts_datasets = []
    if desired!=None:
        for noise_ind in noise_inds:
            TS_y, TS_X = head2head(train_y, train_X, desired=desired, noise=noise_ind)
            ts_datasets.append(torch.utils.data.TensorDataset(TS_X, TS_y))
    else:
        ts_datasets.append(torch.utils.data.TensorDataset(train_X, train_y))

    # Concatenate all datasets into one training set
    training_set = torch.utils.data.ConcatDataset(ts_datasets)
    
    return training_set

#################################################################################################################################################################################

def QuickNN(training_set, n, num_epochs, batch_size, h1_size=0,
            l1_lambda = 0.0, stopEarly = 10, visualize=False):
    # Set the device to use (GPU if available, otherwise CPU)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    # Step 2: Define your neural network
    if h1_size < 1:
        h1_size = n*2
    # Instantiate the neural network model
    inputs = training_set[0][0].shape[0]
    model = Net(input_size=inputs, num_classes=n,
                h1_size=h1_size,
                ).to(device)

    # Define the loss function and optimizer
    criterion = nn.CrossEntropyLoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001, weight_decay=l1_lambda, amsgrad=True)

    # Define the learning rate scheduler
    lr_scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', 
                                                          factor=0.05, patience=5)
    
    # Step 3: Set up your training loop
    train_loader = DataLoader(training_set, batch_size=batch_size, shuffle=True)
    best_loss = float('inf')  # initialize the best validation loss
    early_stop_counter = 0  # initialize the early stopping counter

    if stopEarly > 0:
        print("Early Stopping Initialized")
        # Create the validation set
        val_size = int(len(training_set) * 0.2) # Use 20% of the training set for validation
        val_set, train_set = random_split(training_set, [val_size, len(training_set) - val_size])
        val_loader = DataLoader(val_set, batch_size=batch_size, shuffle=True)


    for epoch in range(num_epochs):
        epoch_loss = 0.0
        for batch_idx, (inputs, targets) in enumerate(train_loader):
            if len(inputs) == 0:
                continue
            inputs, targets = inputs.to(device), targets.to(device)

            # Forward pass
            outputs = model(inputs)
            loss = criterion(outputs, targets)
            epoch_loss += loss.item()

            # Backward Pass
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        if stopEarly > 0:
            # Define the learning rate scheduler
            lr_scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', 
                                                                factor=0.05, patience=5)
            # Evaluate the model on the validation set
            with torch.no_grad():
                val_loss = 0.0
                for inputs, targets in val_loader:
                    inputs, targets = inputs.to(device), targets.to(device)
                    outputs = model(inputs)
                    val_loss += criterion(outputs, targets).item()
            val_loss /= len(val_loader)

            # Check if the validation loss has improved
            if np.round(val_loss,5) < best_loss:
                best_loss = val_loss
                early_stop_counter = 0
            else:
                early_stop_counter += 1
                if early_stop_counter >= stopEarly:  # if the validation loss hasn't improved for 10 epochs, stop training
                    print(f"Early stopping at: Epoch {epoch+1}/{num_epochs}, Loss: {loss.item():.4f}, Learning rate: {optimizer.param_groups[0]['lr']:.6f}")
                    break
        else:
            # Update the learning rate using the scheduler
            lr_scheduler.step(loss)


        if visualize:
            # Print the training loss and learning rate after every epoch
            print(f"Epoch {epoch+1}/{num_epochs}, Loss: {loss.item():.4f}, Learning rate: {optimizer.param_groups[0]['lr']:.6f}")
        
    return model

#################################################################################################################################################################################

def TestModel(test_X, test_y, model, visualize=True):
    # Set the device to use (GPU if available, otherwise CPU)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    test_X = test_X.to(device)
    test_y = test_y.to(device)

    # Evaluate the model on the test set
    with torch.no_grad():
        outputs = model(test_X)
        _, predicted = torch.max(outputs.data, 1)

    results = pd.DataFrame()
    for i in range(min(y), max(y)+1):
        cells = sum(test_y==i).item()
        test_y_i = test_y==i
        y_pred_i = predicted==i
        TP = sum((test_y_i==1) & (y_pred_i==1)).item()
        FP = sum((test_y_i==0) & (y_pred_i==1)).item()
        TN = sum((test_y_i==0) & (y_pred_i==0)).item()
        FN = sum((test_y_i==1) & (y_pred_i==0)).item()
        TPR = TP / np.where(TP+FN == 0, np.nan, TP+FN)
        TNR = TN / np.where(TN+FP == 0, np.nan, TN+FP)
        Prec = TP / np.where(TP+FP == 0, np.nan, TP+FP)
        Accuracy = (TP+TN) / np.where(TP+FP+FN+TN == 0, np.nan, TP+FP+FN+TN)

        res_i = {'Cluster' : le.inverse_transform([i])[0],
            'cells' : cells,
            'TP' : TP,
            'FP' : FP,
            'TN' : TN,
            'FN' : FN,
            'TPR' : TPR,
            'TNR' : TNR,
            'Prec' : Prec,
            'Accuracy' : Accuracy}
        
        res_i = pd.DataFrame([res_i])
        results = pd.concat([results,res_i], ignore_index=True)
        
    if visualize:
        print(results.sort_values(by=['Cluster']))

    return results

#################################################################################################################################################################################
def compute_feature_importance(model, input_data, target_category):
    input_data.requires_grad = True # tell PyTorch to compute gradients with respect to the input
    model.zero_grad()
    output = model(input_data)
    # compute the negative log likelihood loss between the output and the target category
    loss = nn.functional.nll_loss(output, target_category) 
    # compute the gradients of the loss with respect to the input.
    loss.backward()
    # feature importance as the mean absolute value of the gradients over the batch dimension (i.e., over all input examples).
    feature_importance = input_data.grad.abs().mean(dim=0)
    return feature_importance.to('cpu')

def geneSelector(X, genes,device=device):
    X = X.to(device)
    # copy the values from X into X_zeros for the columns we want to keep
    X_sub = X[:, genes]
    return X_sub.to(device)

############################################# CORE CODE ##########################################################################################
##################################################################################################################################################

############################################# LOAD DATA ##########################################################################################
##################################################################################################################################################


# Load data
rdata = pyreadr.read_r(path) # WS


##################################################################################################################################################
# Load And Encode Full Data Set
df = rdata['df']
df['Cluster'] = df['Cluster'].apply(lambda x: x if len(x.split('_')[0]) == 2 else '0' + x) # Standardize cluster names

# Encode the categoric response 
le = LabelEncoder()
df['Cluster'] = le.fit_transform(df['Cluster'])

# Move the response to the end for simply manipulation
cluster_col = df.pop('Cluster')
df.insert(len(df.columns), 'Cluster', cluster_col)

# Remove RBPMS negative cells
df = df[df['Rbpms'] != 0]
df.drop('Rbpms', axis=1, inplace=True)

# Shuffle the data
df = shuffle(df, random_state=42)

# Split the data into input features and labels
X = df.iloc[:, :-1].values.astype(np.float32)
y = df.iloc[:, -1].values.astype(np.compat.long)

# Convert data to PyTorch tensors
X = torch.from_numpy(X)
y = torch.from_numpy(y)

# Split the data into training and test sets
train_size = int(0.8 * len(df))
train_X, test_X = X[:train_size], X[train_size:]
train_y, test_y = y[:train_size], y[train_size:]

# Convert train_X and test_X to PyTorch tensors on the GPU
train_X = train_X.to(device)
test_X = test_X.to(device)

# Compute total number of classes present
n = len(np.unique(y))

##################################################################################################################################################
# Load Tran Hypothesis
hypothesis = rdata['hypothesis']

g1 = hypothesis['Gene1'].copy()
g2 = hypothesis['Gene2'].copy()
g3 = hypothesis['Gene3'].copy()

tran_genes = np.concatenate([g1,g2,g3])
tran_genes = tran_genes[~pd.isna(tran_genes)]

tran_genes = np.where(tran_genes == '4833423E24Rik', 'four833423E24Rik', tran_genes)
tran_genes = np.unique(tran_genes)

tran_gene_indexes = [i for i, col in enumerate(df.columns) if any(gene in col for gene in tran_genes)]

train_X_tran = geneSelector(train_X, tran_gene_indexes)
test_X_tran = geneSelector(test_X, tran_gene_indexes)
##################################################################################################################################################

# Compute Tran Gradients
curriculum_dict['Tran'] = {}
training_set = DesiredGeneCurriculum(train_y, train_X_tran, desired=None)
model = QuickNN(training_set, n, num_epochs, batch_size, stopEarly=0, visualize=False)
results = TestModel(test_X_tran, test_y, model, visualize=False)
curriculum_dict['Tran']['model'] = model
curriculum_dict['Tran']['results'] = results


#################################################################################################################################################################################



# Convert input data to a PyTorch tensor and move it to the GPU
input_data = torch.Tensor(test_X).to(device)


for i in range(0,45):
    curriculum_dict[i] = {}
    training_set = DesiredGeneCurriculum(train_y, train_X, desired=i)
    model = QuickNN(training_set, n, num_epochs, batch_size, stopEarly=5, visualize=False)
    results = TestModel(test_X, test_y, model, visualize=False)
    curriculum_dict[i]['model all genes'] = model
    curriculum_dict[i]['results all genes'] = results
    # Compute feature importance for each target category
    target_category = torch.full((input_data.shape[0],), i, device=device, dtype=torch.long)
    curriculum_dict[i]['feature importance'] = compute_feature_importance(model, input_data, target_category)
    print(f'Curriculum network for cluster {i} complete')

for N in range(2,30):
    # Find the most important features for this number of gradients
    for i in range(0, 45):
        top_n_values, curriculum_dict[i][f'top {N} features'] = torch.topk(curriculum_dict[i]['feature importance'], N, largest=True)
        curriculum_dict[i][f'top {N} genes'] = list(df.columns[curriculum_dict[i][f'top {N} features']])
        print(curriculum_dict[i][f'top {N} genes'])

    top_N_gene_list = [ curriculum_dict[i][f'top {N} genes'] for i in range(0, 45)]

    top_N_gene_list = np.unique(top_N_gene_list)

    print(f' There are {len(top_N_gene_list)} unique genes found by the top {N} gradients from the curriculum models. They are:')
    print(top_N_gene_list)

    common_genes = set(tran_genes).intersection(set(top_N_gene_list))
    print(f'\nThere are {len(common_genes)} genes found by gradient that are in the Tran Proposal they are:')
    x = [print(cg) for cg in common_genes]
    ##################################################################################################################################################

    # Compute Model for the top N genes
    gene_inds_topN = [i for i, col in enumerate(df.columns) if any(g in col for g in top_N_gene_list)]
    train_X_topN = geneSelector(train_X, gene_inds_topN)
    test_X_topN = geneSelector(test_X, gene_inds_topN)

    training_set = DesiredGeneCurriculum(train_y, train_X_topN, desired=None)
    model = QuickNN(training_set, n, num_epochs, batch_size, stopEarly=0, visualize=False)
    results = TestModel(test_X_topN, test_y, model, visualize=True)
    curriculum_dict[f'Model all top {N}'] = {}
    curriculum_dict[f'Model all top {N}']['genes'] = top_N_gene_list
    curriculum_dict[f'Model all top {N}']['model'] = model
    curriculum_dict[f'Model all top {N}']['results'] = results

    ##################################################################################################################################################
    # Compare teh current genes with Tran
    # Recomputed here to prevent in place setting by pandas
    tran_hyp = TestModel(test_X_tran, test_y, curriculum_dict['Tran']['model'], visualize=False)
    gradient = TestModel(test_X_topN, test_y, curriculum_dict[f'Model all top {N}'][f'model'], visualize=False)

    # Add a header indicating which dataframe the row comes from
    gradient.columns = pd.Index([(f'NN Top {N} Gradients', col) for col in gradient.columns])
    tran_hyp.columns = pd.Index([(f'Tran Hypothesis', col) for col in tran_hyp.columns])

    # Concatenate the dataframes horizontally
    results = pd.concat([gradient, tran_hyp], axis=1)

    # Drop the 'Cluster' column from the 'Tran Hypothesis' dataframe
    results = results.drop(('Tran Hypothesis', 'Cluster'), axis=1)

    # Calculate and add the 'TPR Delta' and 'Prec Delta' columns to the 'results' dataframe
    results[('Comparison', 'TPR Delta')] = np.where(results[(f'NN Top {N} Gradients', 'TPR')].isna() | results[('Tran Hypothesis', 'TPR')].isna(),
                                                    np.nan,
                                                    results[(f'NN Top {N} Gradients', 'TPR')] - results[('Tran Hypothesis', 'TPR')])

    results[('Comparison', 'Prec Delta')] = np.where(results[(f'NN Top {N} Gradients', 'Prec')].isna() | results[('Tran Hypothesis', 'Prec')].isna(),
                                                    np.nan,
                                                    results[(f'NN Top {N} Gradients', 'Prec')] - results[('Tran Hypothesis', 'Prec')])

    
    delta_tpr = results[('Comparison', 'TPR Delta')]
    mean_delta_tpr = delta_tpr.mean()

    print(f"Mean TPR Delta for {N} gradients:", mean_delta_tpr)

    delta_prec = results[('Comparison', 'Prec Delta')]
    mean_delta_prec = delta_prec.mean()

    print(f"Mean Prec Delta for {N} gradients:", mean_delta_prec)
    print('#'*100)


##################################################################################################################################################

# Open a file for writing
with open('/home/sam/scRNAseq/RNAscope/RNAscope/curriculum_models_1TO20Gradients.pkl', 'wb') as f:
    # Use pickle to dump the list to the file
    pickle.dump(curriculum_dict, f)