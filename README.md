View under MarkDown mode for better representation.

# Basic Information
* Project: *TDSubgraph*, i.e., *FIDES* in paper
* dev: Renjun Hu 
* homepage: https://hurenjun.github.io/
* Reference: "Fast Computation of Dense Temporal Subgraphs." In *ICDE*, 2017.
* 2017/1/17: reliease the first version of TDSubgraph

# Usage
* **run: compute dense subgraph of temporal graph**
   * E.g.: -run traffic-small.txt 0 10 10 1.0
   * traffic-small.txt: temporal graph
   * 0 & 10: beginning & ending timestamps of temporal graph (from 0 to 10 (including))
   * 10: number k of time intervals that FIDES considers
   * xi: smoothing factor (pg 4 in paper), 1.0 by default

* **runIntv: compute dense subgraph of temporal subgraphs given time intervals (interval.txt)**
   * E.g.: -runIntv traffic-small.txt 0 10 10 1.0
   * Parameter description is the same to *run*
   * 0 & 10 are the beginning and ending timestamps of the temporal graph that we load into main memory. 
   * *It will compute dense subgraph of temporal subgraphs whose beginning and ending timestamps are given in interval.txt.*

* **computeADS: compute dense subgraph on aggregate graph**
   * E.g.: -computeADS traffic-small.txt 0 10
   * Parameter description is the same to *run*
   * *It will compute dense subgraph of aggregate graphs of temporal subgraphs whose beginning and ending timestamps are given in interval.txt*

* **ADSOPT: compute dense subgraph on aggregate graph, with 4 combinations of optimizations**
   * E.g.: -ADSOPT traffic-small.txt 0 10
   * Parameter description is the same to *run*

* **ecrate: compute *ec rate* (i.e., p_{EC}) of temporal graph**
   * E.g.: -ecrate traffic-small.txt 0 10
   * Parameter description is the same to *run*

# Data Format
* temporal graph 
   * each line consists of 4 integers "u,v,t,w"
   * u & v: starting and ending node of the edge
   * t: time stamp
   * w: weight of edge
* interval.txt
   * each line consists of 2 integers "tb,te"
   * tb & te: beginning and ending timestamps of temporal subgraphs

# Baseline & Data
* baseline and synthetic data generator see *meden_1_0_1.zip*
* KPCD groups and ground-truth see *KPCD*
* BJData: available in request
* SynData: generate by synthetic data generator (sizes up to 50G)

For any question, feel free to contact Renjun (hurenjun [AT] buaa.edu.cn)
