import pdb

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

from torchvision import transforms


# %%
class NN(nn.Module):
    def __init__(self, arr=[]):
        super(NN, self).__init__()
        self.relu = nn.ReLU()
        self.fc1 = nn.Linear(30 * 30 * 3, 128)
        self.fc2 = nn.Linear(128, 5)

    def forward(self, x):
        batch_size = x.shape[0]
        x = x.view(batch_size, -1)
        x = self.fc1(x)
        x = self.relu(x)
        x = self.fc2(x)
        return x


# %%
class SimpleCNN(nn.Module):
    def __init__(self, arr=[]):
        super(SimpleCNN, self).__init__()
        self.conv_layer = nn.Conv2d(3, 8, 3)
        self.pool = nn.MaxPool2d(2)
        self.fc1 = nn.Linear(1568, 5)

    def forward(self, x):
        """
        Question 2

        TODO: fill this forward function for data flow
        """
        x = self.conv_layer(x)
        x = F.relu(x)
        x = self.pool(x)
        x = x.view(x.shape[0], -1)
        x = self.fc1(x)
        return x


# %%
basic_transformer = transforms.Compose([transforms.ToTensor()])

"""
Question 3

TODO: Add color normalization to the transformer. For simplicity, let us use 0.5 for mean
      and 0.5 for standard deviation for each color channel.
"""
mean = [0.5, 0.5, 0.5]
standard_deviation = [0.5, 0.5, 0.5]

norm_transformer = transforms.Compose([
    transforms.ToTensor(),
    transforms.Normalize(mean, standard_deviation)])


# %%
class DeepCNN(nn.Module):
    def __init__(self, arr=[]):
        super(DeepCNN, self).__init__()
        """
        Question 4

        TODO: setup the structure of the network
        """

        layers = []
        input_channels = 3

        for i in arr:
            if i == "pool":
                layers.append(nn.MaxPool2d(kernel_size = 2))
            else:
                layers.append(nn.Conv2d(input_channels, i, kernel_size = 3))
                layers.append(nn.ReLU())
            input_channels = i

        self.features = nn.Sequential(*layers)
        self.fc1 = nn.Linear(32 * 12 * 12, 5)

    def forward(self, x):
        """
        Question 4

        TODO: setup the flow of data (tensor)
        """
        x = self.features(x)
        x = x.view(x.shape[0], -1)
        x = self.fc1(x)

        return x


# %%
"""
Question 5

TODO:
    change the aug_transformer to a tranformer with random horizontal flip
    and random affine transformation

    1. It should randomly flip the image horizontally with probability 50%
    2. It should apply random affine transformation to the image, which randomly rotate the image 
        within 5 degrees, and shear the image within 10 degrees.
    3. It should include color normalization after data augmentation. Similar to question 3.
"""

"""Add random data augmentation to the transformer"""
aug_transformer = transforms.Compose([
    transforms.RandomHorizontalFlip(p = 0.5),
    transforms.RandomAffine(degrees = 5, shear = 10),
    transforms.ToTensor(),
    transforms.Normalize(mean, standard_deviation)
])


