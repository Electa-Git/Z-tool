"""
Functions to convert the addmittance matrices from one frame to another
"""
"""
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

__all__ = ['dq_lag2dq_lead','dcdq_lag2dcdq_lead','dq2MSD','dcdq2MSD']

import numpy as np  # Numerical python functions

def dq_lag2dq_lead(Y_old_frame=None, frequencies=None, results_folder=None, file_name=None):
    # Transformation from dq-frame matrices with q-axis lagging to leading
    # The rotation matrix is the product of the abc-to-dq0 transform with q-axis lagging (T1) times the
    # dq0-to-abc transform with q-axis leading (T2, i.e. inverse of the abc-to-dq0 transform with q-axis leading)
    # and then negglecting the zero-sequence component (last row and column): R = T1 @ T2
    R = np.array([[1, 0], [0, -1]]) 
    Y_new_frame = np.matmul(R,np.matmul(Y_old_frame, R))  # = R @ Y @ R is equivalent to multiplying the off-diagonals of Y by -1
    # Save the new matrix
    if file_name is not None and results_folder is not None and frequencies is not None:
        np.savetxt(results_folder+'\\'+file_name+'_dq_leading.txt', np.c_[frequencies,Y_new_frame.reshape(Y_new_frame.shape[0], -1)],
                   delimiter='\t', header="f\t"+"\t".join(["dd","dq","qd","qq"]), comments='')
    return Y_new_frame

def dq2MSD(Y_old_frame=None, frequencies=None, results_folder=None, file_name=None, q_lagging=True):
    # Transformation from a dq-frame matrix to the Modified Sequence Domain (MSD) as described in
    # A. Rygg, M. Molinas, C. Zhang and X. Cai, "A Modified Sequence-Domain Impedance Definition and Its Equivalence to
    # the dq-Domain Impedance Definition for the Stability Analysis of AC Power Electronic Systems," in IEEE Journal of
    # Emerging and Selected Topics in Power Electronics, vol. 4, no. 4, pp. 1383-1396, Dec. 2016, doi: 10.1109/JESTPE.2016.2588733.

    # If the given matrix has q-axis lagging then it is firstly rotated to q-axis leading to match the reference formula
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

def dcdq_lag2dcdq_lead(Y_old_frame=None, frequencies=None, results_folder=None, file_name=None):
    # Transformation from dc,d,q-frame matrices with q-axis lagging to q-axis leading
    # The rotation matrix is the product of the abc-to-dq0 transform with q-axis lagging (T1) times the
    # dq0-to-abc transform with q-axis leading (T2, i.e. inverse of the abc-to-dq0 transform with q-axis leading)
    # then negglecting the zero-sequence component (last row and column), R = T1 @ T2, and considering dc-variables unaffected
    R = np.array([[1, 0, 0],[0, 1, 0], [0, 0, -1]]) 
    Y_new_frame = np.matmul(R,np.matmul(Y_old_frame, R))  # = R @ Y @ R is equivalent to multiplying the q-axis off-diagonals of Y by -1
    # Save the new matrix
    if file_name is not None and results_folder is not None and frequencies is not None:
        np.savetxt(results_folder+'\\'+file_name+'_dq_leading.txt', np.c_[frequencies,Y_new_frame.reshape(Y_new_frame.shape[0], -1)],
                   delimiter='\t', header ="f\tdc-dc\tdc-d\tdc-q\td-dc\td-d\td-q\tq-dc\tq-d\tq-q", comments='')
    return Y_new_frame

def dcdq2MSD(Y_old_frame=None, frequencies=None, results_folder=None, file_name=None, q_lagging=True):
    # Transformation from dc,d,q-frame admittance to the Modified Sequence Domain (MSD) as described in
    # S. Shah and L. Parsa, "Sequence domain transfer matrix model of three-phase voltage source converters," 2016 IEEE
    # Power and Energy Society General Meeting (PESGM), Boston, MA, USA, 2016, pp. 1-5, doi: 10.1109/PESGM.2016.7742009.

    # If the given matrix has q-axis lagging then it is firstly rotated to q-axis leading to match the reference formula
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



dq_lag2dq_lead.__doc__ = """
 Transformation from dq-frame matrices with q-axis lagging to leading
 The rotation matrix is the product of the abc-to-dq0 transform with q-axis lagging (T1) times
 the dq0-to-abc transform with q-axis leading (T2, i.e. inverse of the abc-to-dq0 transform with
 q-axis leading) and then negglecting the zero-sequence component (last row and column): R = T1 @ T2

Required arguments
        Y_old_frame     Numpy matrix of Nx2x2 to be transformed where N is the number of frequency points

Optional arguments
        frequencies     Numpy vector of frequency points. Only needed if the result is to be saved.
        results_folder  Full path of the directory where the result is to be saved.
        file_name       Root name of the saved text file. By default the code adds "dq_leading" to indicate that the matrix has been transformed.
"""

dcdq_lag2dcdq_lead.__doc__ = """
 Transformation from dc,d,q-frame matrices with q-axis lagging to q-axis leading
 The rotation matrix is the product of the abc-to-dq0 transform with q-axis lagging (T1) times
 the dq0-to-abc transform with q-axis leading (T2, i.e. inverse of the abc-to-dq0 transform with
 q-axis leading), then negglecting the zero-sequence component (last row and column): R = T1 @ T2
 while the dc-variables are not afected.

Required arguments
        Y_old_frame     Numpy matrix of Nx3x3 to be transformed where N is the number of frequency points

Optional arguments
        frequencies     Numpy vector of frequency points. Only needed if the result is to be saved.
        results_folder  Full path of the directory where the result is to be saved.
        file_name       Root name of the saved text file. By default the code adds "dq_leading" to indicate that the matrix has been transformed.
"""

dq2MSD.__doc__ = """
Transformation from dq-frame matrix to the Modified Sequence Domain (MSD) as described in
A. Rygg, M. Molinas, C. Zhang and X. Cai, "A Modified Sequence-Domain Impedance Definition and Its Equivalence to
the dq-Domain Impedance Definition for the Stability Analysis of AC Power Electronic Systems," in IEEE Journal of
Emerging and Selected Topics in Power Electronics, vol. 4, no. 4, pp. 1383-1396, Dec. 2016, doi: 10.1109/JESTPE.2016.2588733.
If the given matrix has q-axis lagging then it is firstly rotated to q-axis leading to match the referenced formula

Required arguments
        Y_old_frame     Numpy matrix of Nx2x2 to be transformed where N is the number of frequency points

Optional arguments
        frequencies     Numpy vector of frequency points. Only needed if the result is to be saved.
        results_folder  Full path of the directory where the result is to be saved.
        file_name       Root name of the saved text file. By default the code adds "_MSD_from_dq" to indicate that the matrix has been transformed.
"""

dcdq2MSD.__doc__ = """
Transformation from dc,d,q-frame admittance to the Modified Sequence Domain (MSD) as described in
S. Shah and L. Parsa, "Sequence domain transfer matrix model of three-phase voltage source converters," 2016 IEEE
Power and Energy Society General Meeting (PESGM), Boston, MA, USA, 2016, pp. 1-5, doi: 10.1109/PESGM.2016.7742009.
If the given matrix has q-axis lagging then it is firstly rotated to q-axis leading to match the reference formula

Required arguments
        Y_old_frame     Numpy matrix of Nx2x2 to be transformed where N is the number of frequency points

Optional arguments
        frequencies     Numpy vector of frequency points. Only needed if the result is to be saved.
        results_folder  Full path of the directory where the result is to be saved.
        file_name       Root name of the saved text file. By default the code adds "dq_leading" to indicate that the matrix has been transformed.
"""