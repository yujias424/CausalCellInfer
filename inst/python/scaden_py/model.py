import torch
import random
import warnings
import numpy as np
from tqdm import tqdm
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
warnings.filterwarnings("ignore")

class simdatset(Dataset):
    def __init__(self, X, Y):
        self.X = X
        self.Y = Y

    def __len__(self):
        return len(self.X)

    def __getitem__(self, index):
        x = torch.from_numpy(self.X[index]).float().to(device)
        y = torch.from_numpy(self.Y[index]).float().to(device)
        return x, y

class MLP(nn.Module):
    def __init__(self, input_dim, output_dim, hidden_units, dropout_rates):
        super().__init__()
        self.hidden_units = hidden_units
        self.dropout_rates = dropout_rates
        self.input_dim = input_dim
        self.output_dim = output_dim
        self.model = self._mlp()

    def forward(self, x):
        # x: (n sample, m gene)
        # output: (n sample, k cell proportions)
        return self.model(x)
    # def test(self, x):
    #     # x: (n sample, m gene)
    #     # output: (n sample, k cell proportions)
    #     return self.model(x)

    def _mlp(self):
        mlp = nn.Sequential(nn.Linear(self.input_dim, self.hidden_units[0]),
                            nn.Dropout(self.dropout_rates[0]),
                            nn.ReLU(),
                            nn.Linear(self.hidden_units[0], self.hidden_units[1]),
                            nn.Dropout(self.dropout_rates[1]),
                            nn.ReLU(),
                            nn.Linear(self.hidden_units[1], self.hidden_units[2]),
                            nn.Dropout(self.dropout_rates[2]),
                            nn.ReLU(),
                            nn.Linear(self.hidden_units[2], self.hidden_units[3]),
                            nn.Dropout(self.dropout_rates[3]),
                            nn.ReLU(),
                            nn.Linear(self.hidden_units[3], self.output_dim),
                            nn.Softmax(dim=1)
                            )
        return mlp

def initialize_weight(m):
    if isinstance(m, nn.Linear):
        nn.init.xavier_uniform_(m.weight.data)
        nn.init.constant_(m.bias.data, 0)

class scaden():
    def __init__(self, architectures, train_x, train_y, lr=1e-4, batch_size=128, epochs=20):
        self.architectures = architectures
        self.model512 = None
        self.model256 = None
        self.model1024 = None
        self.lr = lr
        self.batch_size = batch_size
        self.epochs = epochs
        self.inputdim = train_x.shape[1]
        self.outputdim = train_y.shape[1]
        self.train_loader = DataLoader(simdatset(train_x, train_y), batch_size=batch_size, shuffle=True)

    def _subtrain(self, model, optimizer):
        model.train()
        i = 0
        loss = []
        for i in tqdm(range(self.epochs)):
            for data, label in self.train_loader:
                optimizer.zero_grad()
                batch_loss = F.l1_loss(model(data), label)
                batch_loss.backward()
                optimizer.step()
                loss.append(batch_loss.cpu().detach().numpy())
        return model, loss

    def train(self, mode='all'):
        if mode == 'all':
            ##### train
            self.build_model()
            optimizer = torch.optim.Adam(self.model256.parameters(), lr=self.lr, eps=1e-07)
            print('train model256 now')
            self.model256, loss = self._subtrain(self.model256, optimizer)

            optimizer = torch.optim.Adam(self.model512.parameters(), lr=self.lr, eps=1e-07)
            print('train model512 now')
            self.model512, loss = self._subtrain(self.model512, optimizer)

            optimizer = torch.optim.Adam(self.model1024.parameters(), lr=self.lr, eps=1e-07)
            print('train model1024 now')
            self.model1024, loss = self._subtrain(self.model1024, optimizer)

            print('Training of Scaden is done')

    def build_model(self, mode='all'):
        if mode == 'all':
            self.model256 = MLP(self.inputdim, self.outputdim, self.architectures['m256'][0],
                                self.architectures['m256'][1])
            self.model512 = MLP(self.inputdim, self.outputdim, self.architectures['m512'][0],
                                self.architectures['m512'][1])
            self.model1024 = MLP(self.inputdim, self.outputdim, self.architectures['m1024'][0],
                                 self.architectures['m1024'][1])
            self.model1024 = self.model1024.to(device)
            self.model512 = self.model512.to(device)
            self.model256 = self.model256.to(device)
            self.model256.apply(initialize_weight)
            self.model512.apply(initialize_weight)
            self.model1024.apply(initialize_weight)

    def predict(self, test_x, mode='all'):
        test_x = torch.from_numpy(test_x).to(device).float()
        if mode == 'all':
            self.model256.eval()
            self.model512.eval()
            self.model1024.eval()
        if mode == 'all':
            pred = (self.model256(test_x) + self.model512(test_x) + self.model1024(test_x)) / 3
        return pred.cpu().detach().numpy()