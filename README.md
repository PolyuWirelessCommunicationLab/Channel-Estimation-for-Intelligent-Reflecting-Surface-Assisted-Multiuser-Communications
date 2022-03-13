# Channel-Estimation-for-Intelligent-Reflecting-Surface-Assisted-Multiuser-Communications
This code is for paper: [Z. Wang, L. Liu, and S. Cui, "Channel estimation for intelligent reflecting surface assisted multiuser communications: framework, algorithms, and analysis," IEEE Trans. Wireless Commun., vol. 19, no. 10, pp. 6607-6620, Oct. 2020.](https://ieeexplore.ieee.org/abstract/document/9130088)
# Abstract
In intelligent reflecting surface (IRS) assisted communication systems, the acquisition of channel state information (CSI) is a crucial impediment for achieving the beamforming gain of IRS because of the considerable overhead required for channel estimation. Specifically, under the current beamforming design for IRS-assisted communications, $KMN+KM$ channel coefficients should be estimated, where $K$, $N$ and $M$ denote the numbers of users, IRS reflecting elements, and antennas at the base station (BS), respectively. These numbers can be extremely large in practice considering the current trend of massive MIMO (multiple-input multiple-output), i.e., a large $M$, and massive connectivity, i.e., a large $K$. To accurately estimate such a large number of channel coefficients within a short time interval, we propose a novel three-phase pilot-based channel estimation framework in this paper for IRS-assisted uplink multiuser communications, in which the direct channels from the users to the BS, the IRS reflected channels of a typical user, and the IRS reflected channels of the other users are estimated in three consecutive phases, respectively. Under this framework, we analytically prove that a time duration consisting of $K+N+\max(K-1,\lceil (K-1)N/M \rceil)$ pilot symbols is sufficient for the BS to perfectly recover all the $KMN+KM$ channel coefficients for the case without receiver noise at the BS. In contrast to the channel estimation for conventional uplink communications without IRS where the minimum channel estimation time is independent of the number of receive antennas at the BS, our result reveals the crucial role of massive MIMO in reducing the channel estimation time for IRS-assisted communications. Further, for the case with receiver noise, the user pilot sequences, IRS reflecting coefficients, and BS linear minimum mean-squared error (LMMSE) channel estimators are characterized in closed-form, and the corresponding estimation mean-squared error (MSE) is quantified.
# Citation
@article{liu2020IRSchannelestimation,<br> 
  title={Channel estimation for intelligent reflecting surface assisted multiuser communications: framework, algorithms, and analysis},<br> 
  author={Zhaorui Wang and Liang Liu and Shuguang Cui},<br> 
  journal={IEEE Trans. Wireless Commun.},<br> 
  volume={19},<br> 
  number={10},<br> 
  pages={6607--6620},<br> 
  month={Oct.},<br>
  year={2020},<br> 
  publisher={IEEE},<br> 
}
# Note
The code is provided for the benefit of better understanding the results, and is not meant to be used in production.
