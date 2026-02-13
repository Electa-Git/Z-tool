"""
Functions to convert the addmittance matrices from one frame to another
"""
"""
Copyright (C) 2026  Francisco Javier Cifuentes Garcia

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

__all__ = ['dq_lag2dq_lead','dcdq_lag2dcdq_lead','dq2MSD','dcdq2MSD','dq2ab','ab2pn','dq2Jac']

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

def dq2Jac(Y_old_frame=None, frequencies=None, results_folder=None, file_name=None, q_lagging=True, p0=None, q0=None, v0=1.0, v0_dq=None, i0_dq=None):
    # Transformation from a dq-frame admittance to the frequency-domain power flow Jacobian adapted to dq-frame with q-axis lagging from 
    # K. Dey and A. M. Kulkarni, "Passivity-Based Decentralized Criteria for Small-Signal Stability of Power Systems With Converter-Interfaced Generation,"
    # in IEEE Transactions on Power Systems, vol. 38, no. 3, pp. 2820-2833, May 2023, doi: 10.1109/TPWRS.2022.3192302.

    # If the given matrix has q-axis leading then it is firstly rotated to q-axis lagging to match the reference formula
    if not q_lagging:
        R = np.array([[1, 0], [0, -1]])
        Y_old_frame = np.matmul(R,np.matmul(Y_old_frame, R)) 
    
    W = np.array([[0, 1], [-1, 0]])
    if v0_dq is None or i0_dq is None:
        if p0 is None or q0 is None:
            s0 = np.array([[0, 0], [0, 0]])
        else:
            s0 = np.array([[-q0, p0], [p0, q0]])
        v0 = 1 if v0 is None else v0
        E = v0*np.eye(2)
        F = v0*W
        Y_new_frame = np.matmul(np.matmul(E,Y_old_frame), F) + s0
    else:
        E = v0_dq[0]*np.eye(2) + v0_dq[1]*W
        C = np.array([[i0_dq[0], i0_dq[1]], [i0_dq[1], -i0_dq[0]]])
        F = v0_dq[0]*W + v0_dq[1]*np.eye(2)
        Y_new_frame = 3/2*np.matmul(np.matmul(E,Y_old_frame) + C, F)
    
    # Save the new matrix
    if file_name is not None and results_folder is not None and frequencies is not None:
        np.savetxt(results_folder+'\\'+file_name+'_Jac_from_dq.txt', np.c_[frequencies,Y_new_frame.reshape(Y_new_frame.shape[0], -1)],
                   delimiter='\t', header="f\t"+"\t".join(["theta2P","V2P","theta2Q","V2Q"]), comments='')
    return Y_new_frame

def ab2pn(Y_old_frame=None, frequencies=None, results_folder=None, file_name=None, q_lagging=True):
    R = np.eye(2)
    if not q_lagging:
        R = np.array([[1, 0], [0, -1]])   
    A = 0.5 * np.array([[1, 1j], [1, -1j]])  # Same rotation matrix as dq2MSD()
    A_inv = R @ np.array([[1, 1], [-1j, 1j]])
    Y_new_frame = np.matmul(A, np.matmul(Y_old_frame, A_inv))  # Y' = A @ Y @ A^(-1)
    # Save the new matrix
    if file_name is not None and results_folder is not None and frequencies is not None:
        np.savetxt(results_folder + '\\' + file_name + '_pn_from_ab.txt',np.c_[frequencies, Y_new_frame.reshape(Y_new_frame.shape[0], -1)],
                   delimiter='\t', header="f\t" + "\t".join(["pp", "pn", "np", "nn"]), comments='')
    return Y_new_frame

def dq2ab(Y_old_frame=None, frequencies=None, results_folder=None, file_name=None, q_lagging=True, f0=50.0, interpolate=True):
    if np.abs(frequencies[-1] - frequencies[0]) <= f0:
        # Not possible to compute the frame conversion as there are not enough separation between minimum and maximum frequencies
        print(" Matrix frame conversion not possible due to insufficient frequency points")
        return False

    # Transform q-axis from lagging to leading if needed
    if q_lagging:
        R = np.array([[1, 0], [0, -1]])
        Y_old_frame = np.matmul(R, np.matmul(Y_old_frame, R))
    
    # Define the transformation matrices considering q-axis leading
    T = 0.5 * np.array([[1, 1j], [-1j, 1]])
    T_conj = 0.5 * np.array([[1, -1j], [1j, 1]])

    frequencies_new_frame = []
    Y_new_frame = []

    for i, f in enumerate(frequencies):
        if (frequencies[-1] >= np.abs(f-2*f0) >= frequencies[0]):
            if not interpolate:
                if np.abs(f-2*f0) not in frequencies:
                    continue

            frequencies_new_frame.append(f-f0)
            Yplus = Y_old_frame[i]

            idx0 = find_nearest(frequencies, np.abs(f-2*f0))

            Yminus = np.array([
                np.interp(np.abs(f-2*f0), [frequencies[idx0], frequencies[idx0+1]], [Y_old_frame[idx0, 0, 0], Y_old_frame[idx0+1, 0, 0]]),
                np.interp(np.abs(f-2*f0), [frequencies[idx0], frequencies[idx0+1]], [Y_old_frame[idx0, 0, 1], Y_old_frame[idx0+1, 0, 1]]),
                np.interp(np.abs(f-2*f0), [frequencies[idx0], frequencies[idx0+1]], [Y_old_frame[idx0, 1, 0], Y_old_frame[idx0+1, 1, 0]]),
                np.interp(np.abs(f-2*f0), [frequencies[idx0], frequencies[idx0+1]], [Y_old_frame[idx0, 1, 1], Y_old_frame[idx0+1, 1, 1]])
            ])

            # Take the conjugate value for frequencies below f0
            if f-2*f0 < 0:
                for j in range(len(Yminus)): Yminus[j] = Yminus[j].conjugate()
            Yminus = Yminus.reshape(2, 2)
            Y_new_frame.append(np.matmul(T, np.matmul(Yminus, T)) + np.matmul(T_conj, np.matmul(Yplus, T_conj)))
        
        # Estimate the value for the admittance if f-2*f0 is exactly 0
        # elif f-2*f0 == 0 and base_freq:
        #     frequencies_new_frame.append(f-f0)
        #     Yplus = Y_old_frame[i]
        #     idx0 = find_nearest(frequencies, f-2*f0)
        #     if idx0 < 0:
        #         Yminus = Y_old_frame[0] 
        #     else:
        #         Yminus = Y_old_frame[idx0+1]
        #     Y_new_frame.append(np.matmul(T, np.matmul(Yminus, T)) + np.matmul(T_conj, np.matmul(Yplus, T_conj)))

    sorting_idx = np.argsort(frequencies_new_frame)
    frequencies_new_frame = np.array(frequencies_new_frame)[sorting_idx]
    Y_new_frame = np.array(Y_new_frame)[sorting_idx]
    
    if file_name is not None and results_folder is not None:
        np.savetxt(results_folder + '\\' + file_name + '_ab_from_dq.txt',np.c_[frequencies_new_frame, Y_new_frame.reshape(Y_new_frame.shape[0], -1)],
                   delimiter='\t', header="f\t" + "\t".join(["alpha","beta"]), comments='')

    return Y_new_frame, frequencies_new_frame

# Returns index of value nearest to value, but smaller or equal to the original value
def find_nearest(array, value): 
    # Efficient function to find the nearest value to a given one and its position
    idx = np.searchsorted(array, value, side="left")
    if array[idx] == value and idx != len(array)-1:
        return idx
    else:
        return idx-1

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
        q_lagging       Bool to indicate if the q-axis lags the d-axis in the input matrix. Default = True
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
        q_lagging       Bool to indicate if the q-axis lags the d-axis in the input matrix. Default = True
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
        q_lagging       Bool to indicate if the q-axis lags the d-axis in the input matrix. Default = True
        file_name       Root name of the saved text file. By default the code adds "_MSD_from_3x3" to indicate that the matrix has been transformed.
"""

dq2ab.__doc__ = """
Transformation function from d-q frame to the alpha-beta frame.
Required arguments
        Y_old_frame     Numpy matrix of Nx2x2 to be transformed where N is the number of frequency points
        frequencies     Numpy vector of frequency points in the dq frame
        f0              Steady-state frequency in Hz. Default = 50.0

Optional arguments
        results_folder  Full path of the directory where the result is to be saved.
        file_name       Root name of the saved text file. By default the code adds "_ab_from_dq" to indicate that the matrix has been transformed.
        q_lagging       Bool to indicate if the q-axis lags the d-axis in the input matrix. Default = True
        interpolate     Bool to interpolate the matrix to the missing frequencies needed to compute it in the new frame. Default = False.
                        Interpolation leads to more frequency points at the expense of accuracy when the interpolated values are very close to zero.

"""