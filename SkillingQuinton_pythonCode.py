# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 16:16:10 2019

@author: Quinton
"""

###
### Import block
###

#import matplotlib
#matplotlib.use('TkAgg')
#from matplotlib import colors as mcolors
#from matplotlib import pyplot
from math import exp

#import numpy as np
#import pylab as PL
import random as RD
#import scipy as SP
#import networkx as nx
#import time
import sys

###
### Global parameters
###

#connRadius, newSynTimeCutoff (ms), outgoing_halfMax, outgoing_slope, incoming_halfMax, incoming_slope 

stepSize = 0.1

networkSize = 500
targetDensity = 0.1 #This is used to determine grid size
maximalVoltage = 1.1 #Maximum initial voltage
voltageThreshold = 2.
spikeNrns = []
leak = -0.1 #Ionic leakage
connWeight_init = abs(leak) #Initial weight of newly formed connection pulled from a uniform distribution with this as max
learningRate = 0.05 #Strengthening factor during spike-timing-dependent plasticity
newSynTimeCutoff = 2000 #Modulo factor controlling number of steps (technically ms) between growth events
outgoing_halfMax = 5 #Sets half-maximum value of sigmoidal curve used to calculate propensity to form new outgoing connection [# of spikes]
outgoing_slope = outgoing_halfMax/5 #Sets slope of sigmoidal curve used to calculate propensity to form new outgoing connection
incoming_halfMax = 5 #Sets half-maximum value of sigmoidal curve used to calculate propensity to form new incomin connection [# of spikes]
incoming_slope = incoming_halfMax/5 #Sets slope value of sigmoidal curve used to calculate propensity to form new incoming connection
incoming_shape = 0
prevActivity = 0

connRadius = 500#int(sys.argv[2]) #How far away each neuron looks when forming new connections

positionList_x = [] #Global array holding x values of neurons positions
positionList_y = [] #global array holding y values of neuron positions

connectionWeights = {} #Global dictionary holding connection weights, used to easily implement STDP

####
#### Agent and Enviornment Instantiators 
####

class neuron(object):
    """
        This is the class to hold neuron attributes:
        
        Position, activity, etc.
    """
    def _init_(self):
        self.ID = 0
        self.position = []
        self.connections = []
        self.connectionWeights = []
        self.connFormationTimes = []
        self.voltage = RD.random()*voltageThreshold
        self.Input_syn = 0
        self.Input_noise = 0
        self.Input_external = RD.random()*0.8
        self.spikeTimes = []
        self.refractoryTime = -1 #Counter used to track refractory state (>= 0 is in refractory)
        self.prevActivity = 0
        self.neuronsInRange = [] #Tracks the # of neurons in range so as to minimize looping time during connection growth function
        self.data = []
    
def populateSpace_2D():
    """
        This aptly-named function creates a space and populates it with neuron agents
    
    """
    global neurons
    
    grid_size = networkSize/targetDensity
    
    count = 0
    
    for nrn in neurons:
        
        #Generates a new random position
        newPos = [grid_size*RD.random(), grid_size*RD.random()]
        
        nrn.position = newPos
        nrn.ID = count
        
        #Updates global arrays
        positionList_x.append(newPos[0])
        positionList_y.append(newPos[1])
        
        count += 1

####
#### Main Rules of Network Growth
####
        
def outgoingPropensity(act):
    """
        Sigmoidal function to generate probability that a neuron wants to form a new outgoing connection
        
        
        Traditional "S" shape
    """
    
    val = 1./(1. + exp(-(act-outgoing_halfMax)/outgoing_slope))
    return val

def incomingPropensity(act):
    """
        Sigmoidal function to generate probability that a neuron wants to form a new incoming connection
    
        Backward "S" shape
    """
    val = 1.0 - 1./(1. + exp(-(act-incoming_halfMax)/incoming_slope))
    return val
             

def growConnections(isInitial, time_ind):
    
    """
        This function controls the formation of new neurons
        Strengthening of existing neurons happens in a different function
    """
    
    global neurons
    
    grid_size = networkSize/targetDensity
    
    #If "isInitial" is false, new connections are possibly formed
    if(not isInitial):
        
        for nrn in neurons:
            nrn.prevActivity = 0
                    
            #Calculates the activity within the window of length "newSynTimeCutoff" (a parameter)
            for spk in nrn.spikeTimes:
                if((time_ind*stepSize) - spk < newSynTimeCutoff):
                    nrn.prevActivity = nrn.prevActivity + 1
                    
        for nrn in neurons:
              
            outgoingProb = outgoingPropensity(nrn.prevActivity)
            
            possibleConns = []
            
            #Loops through the neurons within range, determined by answer "True" above
            for connNrn in nrn.neuronsInRange:

                #If the connection already exists, it is skipped (assumption of no multi-connections)
                if(connNrn in nrn.connections):
                    continue
                
                #Calculates outgoing and incoming propensities for new connection formation (activity-controlled sigmoidals)
                incomingProb = incomingPropensity(neurons[connNrn].prevActivity)
                        
                randVal = RD.random()
                
                #If a random number on the uniform distribution [0,1] < product of incoming and outgoing propensity, new connection is formed
                if(randVal <= outgoingProb*incomingProb):
                    
                    possibleConns.append(neurons[connNrn].ID)
                    
                    #Updates local connections and strengths
                    #nrn.connections.append(neurons[connNrn].ID)
                    #nrn.connectionWeights.append(RD.random()*connWeight_init)
                    #nrn.connFormationTimes.append(time_ind)
                    
                    #Updates the global strengths which can later be easily modified during STDP
                    #connectionWeights.update({(nrn.ID, neurons[connNrn].ID):nrn.connectionWeights[-1]})
                
            if(len(possibleConns) > 0):
                randConn = int(RD.random()*len(possibleConns))
                nrn.connections.append(neurons[possibleConns[randConn]].ID)
                nrn.connectionWeights.append(RD.random()*connWeight_init)
                nrn.connFormationTimes.append(time_ind)
                    
                #Updates the global strengths which can later be easily modified during STDP
                connectionWeights.update({(nrn.ID, neurons[possibleConns[randConn]].ID):nrn.connectionWeights[-1]})
         
    #If "isInitial" is true, a list of possible neuron targets are generated so as to limit looping times
    else:                
        for nrn in neurons:
            
            nrnPos = nrn.position
            
            #Calculates the maximum possible value within the neurons horizontal connectivity range that it can sense other neurons 
            maxX = nrnPos[0] + connRadius
            if(maxX > grid_size):
                maxX = grid_size
                
            #Same as above but for minimum horizontal distance
            minX = nrnPos[0] - connRadius
            if(minX < 0):
                minX = 0
                
            #Same as above but for maximum vertical distance
            maxY = nrnPos[1] + connRadius
            if(maxY > grid_size):
                maxY = grid_size
            
            #Same as above but for minimum vertical distance
            minY = nrnPos[1] - connRadius
            if(minY < 0):
                minY = 0
                
            #Loop through the global position arrays to determine neurons in range
            for x in range(len(positionList_x)):
                pos_x = positionList_x[x]
                pos_y = positionList_y[x]
                
                #Only the neurons in range can be selected
                if pos_x >= minX and pos_x <= maxX and pos_y >= minY and pos_y <= maxY:
                    
                    #Just a check
                    if(x >= networkSize):
                        break
                    
                    #Assumption: Neurons cannot connect to themselves
                    if(nrn.ID == x):
                        continue
                            
                    #Add the neuron to those in range
                    nrn.neuronsInRange.append(x)
                    
                        
def strengthenConnections(time_ind):
    """
        This function updates the weights of existing connections in an activity-dependent (i.e. spike-timing) way
    """

    global neurons, updateNeurons

    for n_ind in updateNeurons:
        
        nrn = neurons[n_ind]
        
        for connNrn in nrn.connections:
            
            #If the node pointed to by this neuron has not spiked, it is skipped
            #Though connections in the Adjacency matrix A do not necessarily follow Aij = Aji, inactivity is sufficient criteria to skip both Aij and Aji
            if(len(neurons[connNrn].spikeTimes) == 0):
                continue
            
            #Difference in most recent spike times
            spikeTimeDiff = nrn.spikeTimes[-1] - neurons[connNrn].spikeTimes[-1]
            
            #This neuron (nrn) pointing to down-stream target
            preConnTuple = (nrn.ID, neurons[connNrn].ID)
            
            #This neuron (nrn) being pointed to by upstream target
            postConnTuple = (neurons[connNrn].ID, nrn.ID)
            
            #Decreases weight because pre-synaptic neuron spiked after post-synaptic neuron
            if(preConnTuple in connectionWeights):
                
                #Updates the connection strength by decreasing by value given to by exp. decay
                tempVal = connectionWeights[preConnTuple]
                tempVal = tempVal - learningRate*exp(-1.0*spikeTimeDiff/10.0)
                
                if(tempVal < 0.0):
                    tempVal = 0.0
                
                #updates the global connection weight dictionary
                connectionWeights.update({preConnTuple:tempVal})
            
            #Increases weight because post-synaptic neuron spiked after pre-synaptic neuron
            if(postConnTuple in connectionWeights):
                
                #Updates the connection strength by increasing by value given to by exp. decay
                tempVal = connectionWeights[postConnTuple]
                tempVal = tempVal + learningRate*exp(-1.*spikeTimeDiff/10.0)
                
                #Updates the global weight dictionary
                connectionWeights.update({postConnTuple:tempVal})

###
### Main agent state rules
###

def synapticSignal(t_ind, spkTime):
    
    """
        This function generates an action-potential like pulse
    """
    
    timeDiff = stepSize*t_ind - spkTime
    
    #Decay is about 20 ms, same as refractory time
    synOutput = exp(-timeDiff/3.0) - exp(-timeDiff/0.3)
    
    return synOutput
        
def integrateInput(t_ind):
    
    """
        This function calculates the total synaptic input of each neuron
    """
    
    global neurons, updateNeurons
    
    #Loops over neurons that just spiked and sends signal to their downstream targets
    for n_ind in range(len(updateNeurons)):
        
        nrn = neurons[updateNeurons[n_ind]]
        
        #Calculates the signal using the above function
        postSyn_input = synapticSignal(t_ind, nrn.spikeTimes[-1])
        
        for conn_ind in range(len(nrn.connections)):
            
            connNode = neurons[nrn.connections[conn_ind]]
            
            #Updates the synaptic input 
            connNode.Input_syn = connNode.Input_syn + nrn.connectionWeights[conn_ind]*postSyn_input
    
def integrateEquation(t_ind):
    
    """
        This function uses the Euler integration method to update the voltage of each neuron agent
        
    """
    
    global neurons, updateNeurons
    
    for nrn in neurons:
        
        #If a neuron spiked recently, it is in a refractory period where no (realistic) input will cause membrane fluctuations
        if(nrn.refractoryTime >= 0):
            nrn.refractoryTime = nrn.refractoryTime + 1
            if(nrn.refractoryTime == 200):#Fixed refractory time of 20 ms
                nrn.refractoryTime = -1
            
            continue
        
        #Numerical solution to standard Integrate-and-fire differential equation
        #nrn.voltage = nrn.voltage + stepSize*(leak*nrn.voltage + nrn.Input_syn + nrn.Input_noise + nrn.Input_external)
        nrn.voltage = nrn.voltage + stepSize*(leak*nrn.voltage + nrn.Input_external + nrn.Input_syn)
        
        nrn.Input_syn = 0.0
        
        #Check if the neuron spikes and make relevant updates
        if nrn.voltage >= voltageThreshold:
            updateNeurons.append(nrn.ID)
            nrn.voltage = 0.0
            nrn.spikeTimes.append(t_ind*stepSize)
            nrn.refractoryTime = 0
            
    return updateNeurons    

###
### Miscellaneous function
###                   
             
#def draw_adjacency_matrix(G, node_order=None, partitions=[], colors=[]):
#    """
#    - G is a netorkx graph
#    - node_order (optional) is a list of nodes, where each node in G
#          appears exactly once
#    - partitions is a list of node lists, where each node in G appears
#          in exactly one node list
#    - colors is a list of strings indicating what color each
#          partition should be
#    If partitions is specified, the same number of colors needs to be
#    specified.
#    
#    NOT MY FUNCTION
#    TAKEN FROM: http://sociograph.blogspot.com/2012/11/visualizing-adjacency-matrices-in-python.html
#    """
#    adjacency_matrix = nx.to_numpy_matrix(G, dtype=np.bool, nodelist=node_order)
#
#    #Plot adjacency matrix in toned-down black and white
#    fig = pyplot.figure(figsize=(5, 5)) # in inches
#    pyplot.imshow(adjacency_matrix,
#                  cmap="Greys",
#                  interpolation="none")
#    
#    # The rest is just if you have sorted nodes by a partition and want to
#    # highlight the module boundaries
#    assert len(partitions) == len(colors)
#    ax = pyplot.gca()
#    for partition, color in zip(partitions, colors):
#        current_idx = 0
#        for module in partition:
#            ax.add_patch(patches.Rectangle((current_idx, current_idx),
#                                          len(module), # Width
#                                          len(module), # Height
#                                          facecolor="none",
#                                          edgecolor=color,
#                                          linewidth="1"))
#            current_idx += len(module)
#    pyplot.show()
    
def printDataToFiles(fileNum):#, clusterData, clusterData_w):
    
    """
        This function handles printing data to files
    """
    
    global neurons
    
    connFile = open("C:/Users/Quinton/Documents/Classes/ComplexSystem_530/Project/Data/Connectivity_"+str(fileNum)+".dat","w")
    wghtFile = open("C:/Users/Quinton/Documents/Classes/ComplexSystem_530/Project/Data/Weights_"+str(fileNum)+".dat","w")
    timeFile = open("C:/Users/Quinton/Documents/Classes/ComplexSystem_530/Project/Data/connTimes_"+str(fileNum)+".dat","w")
    spikeFile = open("C:/Users/Quinton/Documents/Classes/ComplexSystem_530/Project/Data/spikeData_"+str(fileNum)+".dat","w")
    paramFile = open("C:/Users/Quinton/Documents/Classes/ComplexSystem_530/Project/Data/Parameters_"+str(fileNum)+".dat","w")
    
    #Prints parameters to a parameter file to reduce file name sizes
    #connRadius, newSynTimeCutoff (ms), outgoing_halfMax, outgoing_slope, incoming_halfMax, incoming_slope
    paramFile.write(str(connRadius)+"\n" + str(newSynTimeCutoff)+"\n"+str(outgoing_halfMax)+"\n"+str(outgoing_slope)+"\n"+str(incoming_halfMax)+"\n"+str(incoming_slope))
    
    gCount = 0;
    
    for nrn in neurons:
        
        #Prints the neuron ID, its position, and the unweighted and weighted clustering coefficients to connection and connection weight files
        connFile.write(str(nrn.ID + 1)+"\t" + str(nrn.position[0])+"\t" + str(nrn.position[1]) + "\t" + str(len(nrn.neuronsInRange)))# + "\t" + str(clusterData[gCount]) + "\t" + str(clusterData_w[gCount]))
        wghtFile.write(str(nrn.ID + 1)+"\t" + str(nrn.position[0])+"\t" + str(nrn.position[1]) + "\t" + str(len(nrn.neuronsInRange)))# +  "\t" + str(clusterData[gCount]) + "\t" + str(clusterData_w[gCount]))        
        timeFile.write(str(nrn.ID + 1)+"\t" + str(nrn.position[0])+"\t" + str(nrn.position[1]) + "\t" + str(len(nrn.neuronsInRange)))# +  "\t" + str(clusterData[gCount]) + "\t" + str(clusterData_w[gCount]))
        
        #Determines the maximum size of the loop to reduce the total number of loops
        maxLoopVal = max(len(nrn.spikeTimes), len(nrn.connections))
        
        for ind in range(maxLoopVal):
            
            #If the loop indice is not greater than the number of spikes, those spikes are printed to the spike file
            if(ind < len(nrn.spikeTimes)):
                
                if(nrn.spikeTimes[ind] >= 200):#Cuts out transient
                #Prints neuron ID, its position, and the spike time to the file
                    spikeFile.write(str(nrn.ID + 1)+"\t" + str(nrn.position[0])+"\t" + str(nrn.position[1]) + "\t" + str(nrn.Input_external) + "\t" + str(nrn.spikeTimes[ind]) + "\n")
        
            #Same as above but for the number of connections
            if(ind < len(nrn.connections)):
                #Either prints the post-synaptic neuron ID or the corresponding weight
                connFile.write("\t" + str(nrn.connections[ind] + 1))
                wghtFile.write("\t" + str(connectionWeights[(nrn.ID, neurons[nrn.connections[ind]].ID)]))
                timeFile.write("\t" + str(nrn.connFormationTimes[ind]))
        
        #Prints new lines to the connectivity and weight files
        connFile.write("\n")
        wghtFile.write("\n")
        timeFile.write("\n")
        gCount = gCount + 1
        
    #Closes all the files
    connFile.close()
    wghtFile.close()
    spikeFile.close()
    paramFile.close()
    timeFile.close()
    
def init():
    
    """
        Intialization function
    """
    
    global neurons
    
    positionList_x[:] = []
    positionList_y[:] = []
    connectionWeights.clear()
    
    #Creates neurons
    neurons = []
    for i in range(networkSize):
        neurons.append(neuron())
        neurons[i]._init_()
        neurons[i].ID = i
    
    populateSpace_2D()
    
    #Determine neurons that are in range of one another
    growConnections(True, 0)
    
def runSim(fileNumber):
    
    """
        Main simulation Loop
        
        Progression is:
            1. integrate inputs
            2. Integrate differential equation
            3. 
                a. Grow new connections 
                b. Strengthen existing connections
    """
    
    global neurons, G, updateNeurons
    
    init()
    
    updateNeurons = []
    
    VFile = open("C:/Users/Quinton/Documents/Classes/ComplexSystem_530/Project/Data/Voltage.dat","w")
    
    
    #Time loop
    for time_ind in range(100000):        
    
        integrateInput(time_ind)
    
        #Makes sure that temporary array used to update neurons is empty    
        updateNeurons = []
        
        updateNeurons = integrateEquation(time_ind)
        
        VFile.write(str(time_ind*stepSize) + "\t" + str(neurons[0].voltage) + "\n")
        
        #Assumption: Neurons do not grow connections on every time step and new connections are not formed on the same iteration that existing connections are strengthened
#        if(time_ind > 20000 and time_ind % newSynTimeCutoff == 0):
#            growConnections(False,time_ind)
#        elif(time_ind > 20000):
#            strengthenConnections(time_ind)
#     
        
    #expectedConns = 0    
    #for nrn in neurons:
    #    expectedConns = expectedConns + len(nrn.neuronsInRange)
    
    #print(expectedConns/networkSize)
    #Creates a networkx graph object to easily calculate network properties (clustering coefficient here)
    #G = nx.DiGraph()
    #edgeWeightPairs = []
    
#    for nrn in neurons:
#        G.add_node(nrn.ID)
#        
#        for conn_ind in range(len(nrn.connections)):
#
#            edgeWeightPairs.append((nrn.ID, nrn.connections[conn_ind], nrn.connectionWeights[conn_ind]))
#    
#    #Adds all pre- and post-synaptic neuron IDs and the corresponding pre-to-post connection weight to the graph object
#    G.add_weighted_edges_from(edgeWeightPairs,weight='weight')
#            
#    #Calculates the unweighted and weighted clustering cofficients
#    clusterVals = nx.clustering(G,nodes=None,weight=None)
#    clusterVals_w = nx.clustering(G, nodes=None,weight='weight')
#    
    #Draws the adjacey matrix; not currently being used (done offline)
    #draw_adjacency_matrix(G)
    #PL.show()
    
    #Prints data to files
    printDataToFiles(fileNumber)
    

#Parameter Loops
#connRadius, newSynTimeCutoff (ms), outgoing_halfMax, outgoing_slope, incoming_halfMax, incoming_slope    

fileCount = 1#int(argv[1])

#start = time.time()
#
#for newSynTimeCutoff in [2000]:
#    
#    for outgoing_halfMax in [2.5, 10.5, 18.5, 25.0]:
#        
#        start_b = time.time()
#        
#        for incoming_halfMax in [1, 3, 12, 25.0]:
#            
#            for incoming_shape in [0, 1]:
#            
#                runSim(fileCount)
#                fileCount = fileCount + 1
#        
#        end_b = time.time()
#        
#        print("Loop done in " + str(end_b-start_b) + " seconds...")
#    
#end = time.time()
#print(end-start)

runSim(fileCount)
