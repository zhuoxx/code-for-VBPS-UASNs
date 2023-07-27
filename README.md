# Paper
X. Zhuo, W. Wu, L. Tang, F. Qu, and X. Shen, “Value of information-based packet scheduling scheme for AUV-assisted UASNs,” submitted to IEEE Transactions on Wireless Communications.
# Abstract
In this paper, we propose a value of information (VoI)-based packet scheduling scheme (VBPS) in autonomous underwater vehicle (AUV)-assisted underwater acoustic sensor networks (UASNs), where AUVs act as mobile sensor nodes to collect data from areas not accessible to static nodes and then relay data via static nodes. VoI is a performance metric to measure the importance of data packets with different levels of urgency. The proposed scheme aims to avoid collision with the ongoing packet transmission of static nodes without their accurate global information. In specific, the static node localization stage and the topology construction stage are carried out to obtain the local information. Furthermore, the transmission scheduling stage is implemented to avoid packet collision and formulates a combinatorial optimization problem maximizing VoI under the constraint of packet collision avoidance. To solve this complicated problem, a low-complexity distributed search algorithm is proposed, which exploits the spatial-temporal reuse to establish data packet collision constraints and then determines the next-hop node and data transmission time for AUVs. In addition, a collaborative search algorithm is proposed to avoid packet collision among different AUVs by enabling collaboration among AUVs. Extensive simulation results under various scenarios demonstrate the superior performance of the proposed scheme.
# About this code
Software platform
MATLAB 2019/2020
# Content 
The folder Network performance with different packet arrival rates is provided to produce Fig.5 and Fig.7;  
The folder Network performance with different numbers of AUVs is provided to produce Fig.6;  
The folder Collision probability with vary packet sizes is provided to produce Fig.8.
