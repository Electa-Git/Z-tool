"""
Functions to convert the addmittance matrices from one frame to another

Copyright (C) 2024  Francisco Javier Cifuentes Garcia

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <https://www.gnu.org/licenses/>.
"""

__all__ = ['dq2MSD','dcdq2MSD','ab2pn']

from os import listdir
import numpy as np  # Numerical python functions

def dq2MSD(Y_old_frame=None, frequencies=None, results_folder=None, file_name=None, q_lagging=True):
    # Transformation from dq-frame admittance to the Modified Sequence Domain (MSD) as described in
    # A. Rygg, M. Molinas, C. Zhang and X. Cai, "A Modified Sequence-Domain Impedance Definition and Its Equivalence to
    # the dq-Domain Impedance Definition for the Stability Analysis of AC Power Electronic Systems," in IEEE Journal of
    # Emerging and Selected Topics in Power Electronics, vol. 4, no. 4, pp. 1383-1396, Dec. 2016, doi: 10.1109/JESTPE.2016.2588733.

    # If the given q-axis is lagging then it is firstly rotated to q-axis leading to match the commonly used formulas
    if q_lagging:
        R = np.array([[1, 0], [0, -1]])
    else:
        R = np.eye(2)
    A = 0.5 * np.array([[1, 1j], [1, -1j]]) @ R
    A_inv = R @ np.array([[1, 1], [-1j, 1j]])
    Y_new_frame = np.matmul(A,np.matmul(Y_old_frame, A_inv))  # Y' = A @ Y @ A^(-1)
    # Save the new matrix
    if file_name is not None and results_folder is not None and frequencies is not None:
        np.savetxt(results_folder+'\\'+file_name+'_MSD_from_dq.txt', np.c_[frequencies,Y_new_frame.reshape(Y_new_frame.shape[0], -1)],
                   delimiter='\t', header="f\t"+"\t".join(["pp","pn","np","nn"]), comments='')
    return Y_new_frame

def dcdq2MSD(Y_old_frame=None, frequencies=None, results_folder=None, file_name=None, q_lagging=True):
    # Transformation from dc,d,q-frame admittance to the Modified Sequence Domain (MSD) as described in
    # S. Shah and L. Parsa, "Sequence domain transfer matrix model of three-phase voltage source converters," 2016 IEEE
    # Power and Energy Society General Meeting (PESGM), Boston, MA, USA, 2016, pp. 1-5, doi: 10.1109/PESGM.2016.7742009.

    # If the given q-axis is lagging then it is firstly rotated to q-axis leading to match the commonly used formulas
    if q_lagging:
        R = np.array([[1, 0], [0, -1]])
    else:
        R = np.eye(2)
    A = 0.5 * np.array([[1, 1j], [1, -1j]]) @ R
    A_inv = R @ np.array([[1, 1], [-1j, 1j]])
    A = np.block([[1,np.zeros((1,2))],[np.zeros((2,1)), A]])
    A_inv = np.block([[1,np.zeros((1,2))],[np.zeros((2,1)), A_inv]])
    Y_new_frame = np.matmul(A,np.matmul(Y_old_frame, A_inv))  # Y' = A @ Y @ A^(-1)
    # Save the new matrix
    if file_name is not None and results_folder is not None and frequencies is not None:
        np.savetxt(results_folder+'\\'+file_name+'_MSD_from_3x3.txt', np.c_[frequencies,Y_new_frame.reshape(Y_new_frame.shape[0], -1)],
                   delimiter='\t', header="f\t"+"\t".join(["dc,dc","dc,p","dc,n","p,dc","p,p","p,n","n,dc","n,p","n,n"]), comments='')
    return Y_new_frame

def ab2pn(Y_old_frame=None, frequencies=None, results_folder=None, file_name=None):
    A = 0.5 * np.array([[1, 1j], [1, -1j]])  # Same rotation matrix as dq2MSD()
    A_inv = R @ np.array([[1, 1], [-1j, 1j]])
    Y_new_frame = np.matmul(A, np.matmul(Y_old_frame, A_inv))  # Y' = A @ Y @ A^(-1)
    # Save the new matrix
    if file_name is not None and results_folder is not None and frequencies is not None:
        np.savetxt(results_folder + '\\' + file_name + '_pn_from_ab.txt',np.c_[frequencies, Y_new_frame.reshape(Y_new_frame.shape[0], -1)],
                   delimiter='\t', header="f\t" + "\t".join(["pp", "pn", "np", "nn"]), comments='')
    return Y_new_frame
