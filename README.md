# MSc Thesis: Data Analysis of Particle Collisions at Large Hardon Collider

## Project: Exploring Beyond the Standard Model with High-Energy Physics

### Introduction and Motivation
- High-energy physics endeavours to understand fundamental nature through particle collisions at extreme energies.
- Higgs boson's prediction in the 1960s addressed shortcomings in the Standard Model (SM).
- LHC experiments (ATLAS, CMS) amassed vast data; SM explains only 4% of the universe.
- Pursuit of physics models beyond the SM seeks to unravel the unknown phenomena.

### Project Goal
- Investigate precision in the SM via the LHC's potential.
- Discover novel physics and maintain background rejection.
- Employ machine learning (ML) algorithms to enhance signal efficiency.
- Utilize data simulation tools (MadGraph, Pythia, Delphes) to generate particle collision events.
- Machine learning tools like neural networks and boosted decision trees (BDT) play a crucial role.

## Simulation Process Summary
- Specify collision process: particles, energy, initial/final states.
- MadGraph employs Feynman diagrams for cross sections, generating events with user-defined constraints.
- Pythia simulates parton showering and hadronization after hard scattering.
- Delphes simulates detector response to collisions, handling particle interactions.
- ROOT extracts data, and reveals particle behaviour, and interactions.

## Study of a very specific process for the production of light particles and photons

### Introduction

The pursuit of unravelling the mysteries of the universe has led physicists to embark on groundbreaking experiments at the Large Hadron Collider (LHC). Among the multitude of processes studied, a specific focus lies on the production of light particles and photons. This endeavour delves into the realm of particle interactions, shedding light on the intricate mechanisms that govern our universe. In this study, we delve into the details of this particular process, exploring data generation, transformation, and classification.

### Data Generation and Process Details

The core of our investigation involves simulating events that encapsulate the production of particles and photons. Two pivotal types of events are generated, each holding distinct significance:
- **Background Events**: These events involve the interaction of quarks and antiquarks (q + ¯q) resulting in the production of jets and a photon (γ). The final state comprises two jets and one photon.
- **Signal Events**: In this scenario, quarks and antiquarks again interact, leading to the production of a Z boson (Z) which further decays into a bottom quark (b), an anti-bottom quark (¯b), and a photon. This manifests as bottom, anti-bottom jets, and a photon.

![image](https://github.com/shirsh10mall/MSc-Project---Studies-at-Large-Hadron-Collider/assets/87264071/df0ecc6f-8d1c-4dc2-b9fd-8607c6a857a0)


MadGraph software takes the lead in generating these events, while Pythia software steps in to simulate the subsequent processes of showering and hadronization. The resultant data is then processed for analysis, with a focus on seven critical physical quantities. These include the invariant mass of parent particles of jets, Transverse Momentum (pT) of leading and second-leading jets, energy of these jets, and the energy and mass of radiated photons.

### Transforming Data into Images

To harness the power of deep learning networks, we take a novel approach of transforming our data into images. The Matplotlib library in Python proves instrumental in this endeavor. The dataset, enriched with information about energy, angular positions (△ϕ), η (pseudorapidity), and energy of each constituent within a jet, takes a new form as images. The process involves constructing a grid where η spans one axis, while △ϕ spans the other with distinct upper and lower bounds. Energy summation for constituents is discretized between 0 and 255, effectively assigning color channels. This transformation results in images that capture the essence of particle interactions.

![image](https://github.com/shirsh10mall/MSc-Project---Studies-at-Large-Hadron-Collider/assets/87264071/9cbcf1d4-c581-4b30-b3b8-c7f15a1f7da0)

![image](https://github.com/shirsh10mall/MSc-Project---Studies-at-Large-Hadron-Collider/assets/87264071/67de5ce9-2a07-4466-b2f1-d47fd539a1d9)

![image](https://github.com/shirsh10mall/MSc-Project---Studies-at-Large-Hadron-Collider/assets/87264071/50427a4f-750f-4a3f-9278-2ae406d6f9cc)


![image](https://github.com/shirsh10mall/MSc-Project---Studies-at-Large-Hadron-Collider/assets/87264071/e2d83ed0-f963-4350-8e75-01b120c3b8a9)


## Electron/Photon Classification

The crux of our study revolves around the classification of particles, specifically electrons and photons. The underlying challenge is to discern between signal (electron) and background (photon) processes, crucial for a comprehensive understanding of particle behavior. The dataset encompasses two essential columns:
- **"X" Column**: This holds tensor data with dimensions (249000, 32, 32, 2), encapsulating events and image dimensions. The first two dimensions represent the number of events and the image size (32x32), while the third dimension introduces channels that represent hit energy and time.
- **"y" Column**: This column contains labels, assigning 0 to photons and 1 to electrons. These labels establish two distinct classes that are the focus of our classification task.

![image](https://github.com/shirsh10mall/MSc-Project---Studies-at-Large-Hadron-Collider/assets/87264071/a1292c82-05c1-46f3-85d8-23d12f879109)


In our pursuit to classify these particles, Convolutional Neural Networks (CNNs) emerge as a robust solution. Leveraging the TensorFlow library, we construct a CNN architecture comprising various layers. This intricate architecture is tailored to analyze image data, uncover patterns, and discern the subtle differences between electrons and photons.

In essence, this study encapsulates the journey of exploring particle production, data transformation, and classification at the LHC. By delving into the minutiae of this process, we contribute to the collective understanding of the fundamental building blocks of our universe.

## Regression Model for Transverse Momentum Prediction

### Data Generation and Extraction
For each process, we generated approximately 10,000 events and extracted essential features, including invariant mass, pseudorapidity, and transverse momentum (pT) of the first leading electrons, along with their x and y components (px, py).

### Transverse Momentum (pT) Definition
Transverse momentum (pT) represents the particle's traversal in the z direction within the xz-plane. The transverse plane is the xy-plane, where particle scattering occurs. In the center-of-mass frame, the angle of scattering is θcm, and the azimuthal angle is ϕ. The pT of a particle is defined as pT = pcos(θcm) = px, while the longitudinal momentum is pT = psin(θcm) = pz.

### Training Approach
Our objective was a regression task targeting transverse momentum (ptl1). We utilized both background and signal event datasets. The background event dataset involves p+p → e+ + e− processes, and the signal event dataset involves p+p → Z → e+ + e− processes. Notably, features included px, py, and pT momenta of the first leading electron.

To train our model, we employed an Artificial Neural Network (ANN) on the background event data. Subsequently, this trained model was used for predicting transverse momentum for the signal events. The background dataset was divided into training (3026 examples) and validation (757 examples) subsets, while the signal dataset comprised 4231 examples.

The Neural Network architecture consisted of an input layer, a hidden layer with 100 nodes utilizing a linear activation function, and an output layer with a single node and linear activation. We initialized weights with small random values, and the loss function was set as the mean square error. Our optimizer of choice was the Adam Optimizer.

### Results
Our focus on Root Mean Squared Error (RSME) values yielded promising results. As the number of epochs increased, the RSME for both training and validation data consistently decreased. These RSME values remained low and similar for both sets, indicating effective model tuning and successful avoidance of underfitting or overfitting issues.


## Classification Model for Spin States Prediction

### Data Details
The γγ final state has gained significance post the discovery of the Higgs boson in the di-photon channel. This resonant di-photon channel can be described by the formula: `p + p → χ + x → γγ + x`, where χ is an electrically neutral resonance, either spin 0 or spin 2, and x represents residual protons after spectator quarks are obtained [4]. The presence of x along the line of particles renders them unobserved. While heavy resonance χ may degenerate into nn or nγ, fundamental particle quarks and gluons might degenerate into photons, leading to a state with multiple photons [4]. Discriminating between spins 0 and 2 of χ based on the distributions of Delta-eta (Δη) has been studied [4]. However, this project aims to employ advanced deep learning algorithms for the same purpose.

### Process Equation
The process equation for this study is: `p + p → χspin=0,2 → nn → γγ`. It is challenging to distinguish between spins 0 and 2 of χ due to their identical decay processes. The differences in predicted Δη distributions are used to determine the minimum number of processes required to differentiate between collision events corresponding to various spins. The distribution of Δη for spin 2 states is non-central, while that of spin 0 states is central [4].

### Model Construction
An Artificial Neural Network (ANN) is employed for classifying the two spin states. Data is cleaned and labeled such that Δη values corresponding to spin 2 states are labeled as "1" and those for spin 0 states are labeled as "0". The merged data is split into training and testing datasets, with approximately 75% used for training. The ANN model consists of a single hidden layer with 10 units and ReLU activation. The output layer employs sigmoid activation for classification. The model is trained using the mean squared error loss function and the Adam optimizer.

## Model Performance Assessment

The performance of our model is meticulously evaluated using comprehensive metrics derived from both training and validation sets. These metrics are indicative of the model's efficiency in learning and generalizing the patterns:

**Training Metrics**:
- Loss: 0.5469
- Accuracy: 0.7295
- Area Under Curve (AUC): 0.7964

**Validation Metrics**:
- Loss: 0.5533
- Accuracy: 0.7240
- AUC: 0.7935

### Results
The training and validation loss curves demonstrate a monotonically decreasing trend with increasing epochs, indicating that the model is neither underfitting nor overfitting. The ROC score approaches 1 with increasing epochs, signifying improved classification performance. The AUC score for the training dataset is 0.985, indicating strong performance. However, the AUC score for the test dataset is 0.857, implying a need for model improvement or alternative methods for further enhancement. The ROC curve of the test dataset resembles the random line, suggesting scope for refining the model or adopting different approaches for building a new model.


 -- Simulating proton-proton collision events, generating data of physical quantities using Linux-based simulation software, and constructing new features and images using theoretical learning.
 -- Building Regression predictive models for physical quantities, Classification and Anomaly detection models to classify signal and background events using ANN, CNN and Auto-Encoders. Creating Clustering models for Jets clustering using the Anti-KT Algorithm.

## Project Files: https://csciitd-my.sharepoint.com/:f:/g/personal/phs217221_iitd_ac_in/EoYKZMJhGwpApxWJyp0KYG8BqkFBGLgCjajxn9ivlxKZpg?e=c0h9KZ
