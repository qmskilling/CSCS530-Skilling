# Model Proposal for Investigations of Spatial Constraints on Neuronal Networks and Activity

_Quinton Skilling_

* Course ID: CMPLXSYS 530
* Course Title: COmputer Modeling of Complex Systems
* Term: Winter, 2019

&nbsp;

### Goal
*****

Experimental and computational studies of neural structures and activity abound but many consider neural systems in fully formed networks. This project aims to understand how spatial dimension constrains how neural networks are formed and how those networks shape corresponding neural activity (i.e. "spikes"). This may give rise to a better understanding of how networks form in vivo and in vitro and how that translates to cognition.

&nbsp;
### Justification
******

Neuronal networks are complex, dynamical systems whose behavior emerges from interactions between individual neurons. Neuron shape and function come in different varieties with heterogeneities in the type and number of synaptic connections that form networks. Using agent-based modeling to simulate individual neurons and their connections will allow me to examine how these heterogeneities affect the growth and maintenence of networks and how those networks consequently shape neuronal activity.

&nbsp;
### Main Micro-level Processes and Macro-level Dynamics of Interest
*****

As I will outline below, I want to examine how different spatial dimensions affect the network structure and subsequent activity of those networks. The micro-level processes I will focus on are the rules for synapse formation and the activity of individual neurons. Synapses will form probabilistically by surveying their local environment for (possibly noisy) activity. The activity itself will follow a simple integrate-and-fire formalism, where persistent, environment-driven input will drive early activity and emergent synaptic input will come into effect as networks form. These networks will be probed for structure and other network features including clustering, centrality, and connectivity distribution. Further, by modifying the integrate-and-fire rules to allow for different excitability types, I will examine emergent behavior of network activity by analyzing spiking dynamics for phase coherence/synchronization and effective network connectivity, i.e. networks calculated based on activity instead of structure.

&nbsp;

## Model Outline
*****
&nbsp;

### 1) Environment
The environment of this model will be mostly passive but will be used in conjunction with changing agent parameters to examine properties of structural and effective neural networks.

1. Enviornment will be a grid-free square space
  - Neurons will be placed randomly on the grid
2. Boundary conditions will be parameterized to be wrapped or un-wrapped
  - Initially, only unwrapped boundary conditions will be used.
  - Wrapped boundary conditions will be implemented if time permits.
3. Dimensionality will be a parameter
  - 2D vs 3D effects on network properties will be examined
4. Enivornment-owned variables
  - Cell-packing density
  - Dimensionality
  - Background noise
5. Enviornment-owned methods/procedures
  - Neuron placement
  - Noise <b> will only be set by the enviornment </b> and will be accounted for by an agent-owned procedure.
&nbsp;

### 2) Agents

Agents in this system will the neurons and connections will be sub-agents belonging to neurons.

1. Agent-owned variables
  - Voltage
  - Survey radius and timescale (to form connections)
  - Connectivity (once connections are formed)
  - Intrinsic frequency/firing threshold
  - Excitability type
  
2. Agent-owned methods/procedures
  - Calculate input from connections
  - Calculate input from environment (noise)
  - Update voltage
  - Spike and reset
  - Adjust firing frequency

&nbsp;

### 3) Action and Interaction


#### Interaction Topology

Interactions will be done in three ways. 

1. Neurons will survey their local environment (controlled by <b> survey radius </b>) and will form connections probabilistically, taking into account time-of-activity correlations with nearby neurons.

2. Once connections are formed, neurons will receive input from their pre-synaptic partners based on the connectivity strength between partners.

3. Connectivity strength will be adjusted following a spike-timing-dependent plasticity rule, where pre-before-post firing will increase the strength and post-before-pre will decrease the strength.

#### Action Sequence

Neurons will update synchronously with the sequence on each turn as follows:

1. Neurons calculate their input
2. Neurons integrate their input and potentially activate
3. Connection strengths are adjusted (if they exist)
4. Neurons survey for synapse creation <b> if they did not fire </b>

&nbsp;

### 4) Model Parameters and Initialization

Global parameters include: the spatial dimension of the model (2D vs 3D); the growth rate of connection formation; the growth-rate of connection strengthening/weakening; level of background noise; and the heterogeneity factors in neurons' firing rate/threshold, survey radius and timescale, and excitability type (integrator or resonator).

The model will be initialized by first populating the environment with neurons and then assigning random values for a neuron's voltage and parameters controlled by heterogenous factors (those listed in the last point, above). 

The model procedure is as follows:

1. One-time initialization including steps listed above.

2. Neurons calculate their input
-Synaptic input is calculated by multiplying synaptic weight by an output profile meant to mimick the shape of action potential generation in biophysical neurons
-Noise from the environment is calculated. Input will be given based if the noise of the system is greater than a real-valued number pulled from a uniform distribution between 0 and 1. 
  
3. Neurons integrate their input and potentially activate
-Neurons update their <b> voltage </b> following the Integrate-and-Fire formalism: V(t+1) = V(t) + I_{synaptic} + I_{noise} - I_{leak}
-If <b> V(t) > V(threshold) </b>:
  - V(t) fires and the timing is tabulated
  - V(t) = 0
-If the neuron fired, it is flagged to update synaptic output on the next time step (Step 1)

4. Connection strengths are adjusted (if they exist)
  - Spike-timing-dependent plasticity
    - Pre-before-post will increase connectivity strength
    - Post-before-pre will decrease connectivity strength
5. Neurons survey for synapse creation <b> if they did not fire </b>
  - Neurons correlate their activity to those in <b> survey radius </b>
  - If correlation is high over <b> survey timescale </b>, a connection forms <b> to the surveying neuron

### 5) Assessment and Outcome Measures

Successful implementation of this model will yield quantifiable network structures (actual connections) and population spiking data. 

Network structures will be probed for measures such as component size (is the network fully connected), clustering (number of pairwise connections vs triplet connections), and connectivity distribution, among others; Newman will be a source. 

Spiking data will be analyzed for such things as synchronization and phase coherence, effective network connectivity (e.g. using cross-correlation or mutual information on all pairwise spike timings)

These data will be combined and used to classify the system in terms of the spatial dimension constraint. The idea is to get an idea of what kinds of networks form and what the corresponding activity looks like in terms of different growth environments.
&nbsp;

### 6) Parameter Sweep
