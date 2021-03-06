# Model Proposal for Investigations of Spatial Constraints on Neuronal Networks and Activity

_Quinton Skilling_

* Course ID: CMPLXSYS 530
* Course Title: COmputer Modeling of Complex Systems
* Term: Winter, 2019

&nbsp;

### Goal
*****

Experimental and computational studies of neural structures and activity abound but many consider neural systems in fully formed networks. This project aims to understand how spatial dimension constrains how neural networks are formed and how those networks shape corresponding neural activity (i.e. "spikes"). This may give rise to a better understanding of how networks form in vivo and in vitro and how that translates to neural processes such as memory consolidation.

&nbsp;
### Justification
******

Neuronal networks are complex, dynamical systems whose behavior emerges from interactions between individual neurons. Neuron shape and function come in different varieties with heterogeneities in the type and number of synaptic connections that form networks. Using agent-based modeling to simulate individual neurons and their connections will allow me to examine how these heterogeneities affect the growth and maintenence of networks and how those networks consequently shape neuronal activity.

&nbsp;

*__LS COMMENTS__: Great description and overview of why ABM is being used.*

### Main Micro-level Processes and Macro-level Dynamics of Interest
*****

As I will outline below, I want to examine how different spatial dimensions affect the network structure and subsequent activity of those networks. The micro-level processes I will focus on are the rules for synapse formation and the activity of individual neurons. Synapses will form probabilistically by surveying their local environment for (possibly noisy) activity. The activity itself will follow a simple integrate-and-fire formalism, where persistent, environment-driven input will drive early activity and emergent synaptic input will come into effect as networks form. These networks will be probed for structure and other network features including clustering, centrality, and connectivity distribution. Further, by modifying the integrate-and-fire rules to allow for different excitability types, I will examine emergent behavior of network activity by analyzing spiking dynamics for phase coherence/synchronization and effective network connectivity, i.e. networks calculated based on activity instead of structure.

&nbsp;

*__LS COMMENTS:__ Great explanation.*

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
  - Noise will only be _set_ by the enviornment but will be accounted for by an agent-owned procedure (_Calculate input from environment_).
&nbsp;


Example Code & psuedocode: 

```python
networkSize = 1000
targetDensity = 0.1
backgroundNoise = 0

def populateSpace_2D():

    global neurons
    
    grid_size = networkSize/targetDensity
    
    for nrn in neurons:
        
        newPos = [grid_size*RD.random(), grid_size*RD.random()]
        
        nrn.position = newPos
        
def init():
    global neurons
    
    neurons = []
    for i in range(networkSize):
        neurons.append(neuron(i))
    
    populateSpace_2D()

```

*__LS COMMENTS__: Very solid.*

### 2) Agents

Agents in this system will be the neurons and their connections will be sub-agents.

1. Agent-owned variables
  - Voltage
  - Survey radius and timescale (to form connections)
  - Connectivity strength (once connections are formed)
  - Intrinsic frequency/firing threshold
  - Excitability type
  
2. Agent-owned methods/procedures
  - Calculate input from connections
  - Calculate input from environment (noise)
  - Update voltage
  - Spike and reset
  - Adjust firing frequency
  
  
Example code and psuedocode:

```python
def calculateInput(time_ind):
    
    global neurons
    
    #Loops through neurons
    for nrn in neurons:
        
        #Copies conneciton list
        conns = nrn.connections
        
        #if there are no connections, skip this neuron
        if conns is null:
            continue
        
        #Loops through the connections
        for conn in conns:
            #creates a temp variable for pre spike time
            preSynSpikeTime = neurons[conn].spikeTime.end()
            
            #updates the synaptic input
            #Integrate-and-fire neurons are Type 1 by default, so I will have to think of a succint way to introduce type 2
            nrn.Input_syn += weights[preSyn]*(exp(-(time_ind*stepSize - preSynSpikeTime)/30) - exp(-(time_ind*stepSize - preSynSpikeTime)/3.))

        #updates the noise input based on probability
        if RD.random() < noiseProb:
            nrn.Input_noise += noiseVal
            
def updateVoltage(time_ind):
    
    global neurons
    
    #loops through the neurons
    for nrn in neurons:
    
        #Updates voltage using synaptic input and noise input
        nrn.voltage += stepSize*(nrn.Input_syn + nrn.Input_noise)
        
        count = 1
        
        #spiking criteria: if neuron spikes, tabulate the time, update the spike list, and reset voltage
        if nrn.voltage >= voltageThreshold:

            nrn.voltage = 0
            nrn.spikeTimes.append(time_ind*stepSize)
            spikeNrns.append(count)
            
            count += 1
            

```

&nbsp;

*__LS COMMENTS__: Good.*


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
4. Neurons survey for synapse creation _if they did not fire_

&nbsp;

*__LS COMMENTS__: Great.*

### 4) Model Parameters and Initialization

Global parameters include: the <b> spatial dimension </b> of the model (2D vs 3D); the new <b> connection growth-rate </b>; the <b> connection strength growth-rate </b>; level of background <b> noise </b>; and the heterogeneity factors in neurons' firing rate/threshold, <b> survey radius </b>  and <b> survey timescale </b> , and <b> excitability type </b> (integrator or resonator). Ideally, <b> excitability type </b> can rely on the <b> spatial dimension </b>, mimicking early network formation in vivo (3D, where molecular gradients control position of excitability type) or network formation in vitro (2D, where excitability type is mixed); this may be too complicated, however. 

The model will be initialized by first populating the environment with neurons and then assigning random values for a neuron's <b> voltage </b> and parameters controlled by heterogenous factors (those listed in the last point, above). 

The model procedure is as follows:

1. One-time initialization including steps listed above.

2. Neurons calculate their input
-Synaptic input is calculated by multiplying synaptic weight by an output profile meant to mimick the shape of action potential generation in biophysical neurons
-Noise from the environment is calculated. Input will be given based if the noise of the system is greater than a real-valued number pulled from a uniform distribution between 0 and 1. 
  
3. Neurons integrate their input and potentially activate
  - Neurons update their <b> voltage </b>
    - V(t) = V(t-1) + I_{synaptic} + I_{noise}
  - If V(t) > V(threshold), the neurons fires an action potential (sent to post-synaptic partner, updated in Step 1 on the next tick) and is reset to zero.

4. Connection strengths are adjusted (if they exist)
  - Spike-timing-dependent plasticity
    - Change in spike timing is calculated for each neuron and its synaptic partners,  dt = t(pre spike) - t(post spike)
    - If dt < 0 (pre-synaptic firing before post-synaptic), w(pre-post) increases
    - If dt > 0 (post-synaptic firing before pre-synaptic), w(pre-post) decreases
    - Timings should work that with tonic firing, weights will not change on average
5. Neurons survey for synapse creation _if they did not fire_
  - Neurons correlate their activity to those in <b> survey radius </b> and <b> survey timescale </b>
  - If correlation is high, a connection forms from the _surveyed_ neuron to the _surveying_ neuron

&nbsp;


*__LS COMMENTS__: Extremely clear and well-thought out. Good idea on handling new synapse creation.*

### 5) Assessment and Outcome Measures

Successful implementation of this model will yield quantifiable network structures (actual connections) and population spiking data. 

Network structures will be probed for measures such as component size (is the network fully connected), clustering (number of pairwise connections vs triplet connections), and connectivity distribution, among others; Newman will be a source. 

Spiking data will be analyzed for such things as synchronization, phase coherence, bursting, and effective network connectivity (e.g. using cross-correlation or mutual information on all pairwise spike timings), among others that may arise from the data.

These data will be combined and used to classify the system in terms of the spatial dimension constraint. The idea is to determine of what kinds of networks form and what the corresponding activity looks like in terms of different growth environments and to relate this information back to experimental findings.

&nbsp;


### 6) Parameter Sweep

I plan to sweep through the following parameters:
- <b> survey radius </b>
- <b> survey timescale </b>
- <b> cell-packing density </b>
- <b> background noise </b>
- <b> maximal connection strength </b>
- <b> connection growth rate </b>
- <b> connection strength growth rate </b>

I plan to normalize as many of these parameters as I can. The <b> survey radius </b>, for example, will be normalized by the expected average spatial distance between neurons (based on <b> cell-packing density </b>) and the <b> connection strength growth rate </b> will be normalized by the <b> maximal connection strength </b>.

If normalization is successful, I expect that values close to unity in <b> survey radius </b> and <b> connection strength growth rate </b> will result in fully connected networks and fully synchronous dynamics (in the case of Type 2 excitability), respectively. However, with high levels of noise, I expect disruption to both a fully-connect network structure and asynchronous dynamics. 

I also expect that large values of the <b> survey radius </b> and long <b> survey timescales </b> will increase network connectivity and will serve as a high-frequency perturbation to neural activity.

Overall, I expect large values of the above-listed parameters will result in densely connected networks with activity following the mean population excitability type (i.e., type 2 networks will be synchronous whereas type 1 will not). However, speculating about how intermediate values of these parameters will affect network formation and emergent activity is difficult; this is what I want to characterize with the model.


*__LS COMMENTS__: Fantastic project that has obviously been well-thought out! Very much look forward to seeing what this produces. 19.5/20*
