#!/usr/bin/env python3
"""
Convert 3D velocity field arrays from HDF5 to binary format
compatible with the Phantom and sphNG smoothed particle 
hydrodynamics codes
"""

import h5py
import numpy as np
import struct
import argparse

def write_fortran_unformatted(filename, data):
    """
    Write data in Fortran unformatted format with proper record markers.
    
    Parameters:
    filename: str - output filename
    data: numpy array - data to write
    """
    
    # Convert to 32-bit float (Fortran real(kind=4))
    data_f32 = data.astype(np.float32)
    
    # The Fortran code writes: (((vel(i,j,k,component),i=1,N),j=1,N),k=1,N)
    # This means i varies fastest, then j, then k
    # To match Fortran output exactly, we need to transpose our Python data
    # from [i,j,k] to [k,j,i] ordering before writing
    # Transpose to match Fortran coordinate system: (i,j,k) -> (k,j,i)
    data_transposed = np.transpose(data_f32, (2, 1, 0))

    # Ensure Fortran ordering (column-major, where first index varies fastest)
    data_fortran = np.asfortranarray(data_transposed)

    # Calculate record length in bytes
    record_length = data_fortran.nbytes

    with open(filename, 'wb') as f:
        # Write opening record marker (4-byte integer)
        f.write(struct.pack('i', record_length))

        # Write the data
        f.write(data_fortran.tobytes())

        # Write closing record marker (4-byte integer)
        f.write(struct.pack('i', record_length))

def convert_velocity_field(input_file, output_files):
    """
    Convert HDF5 velocity field to binary files readable by Fortran.
    
    Parameters:
    input_file: str - path to HDF5 file containing vx, vy, vz datasets
    output_files: list - [cube_v1.dat, cube_v2.dat, cube_v3.dat] output file names
    """
    
    # Read HDF5 file
    print(f"Reading {input_file}...")
    with h5py.File(input_file, 'r') as h5file:
        # Read the velocity components
        vx = h5file['vx'][:]
        vy = h5file['vy'][:]
        vz = h5file['vz'][:]
        
        print(f"Array shape: {vx.shape}")
        print(f"Data type: {vx.dtype}")
    
    # List of velocity arrays and corresponding output files
    velocity_arrays = [vx, vy, vz]
    
    for i, (velocity, output_file) in enumerate(zip(velocity_arrays, output_files)):
        print(f"Converting component {i+1} to {output_file}...")
        
        # Write as Fortran unformatted binary file with record markers
        write_fortran_unformatted(output_file, velocity)

def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Convert HDF5 velocity field to Fortran binary files")
    parser.add_argument('input_file',
                        nargs='?',
                        default='velocity_field_seed0.h5',
                        help='HDF5 input file (default: velocity_field_seed0.h5)')

    args = parser.parse_args()

    # Output file names (always the same)
    output_files = ["cube_v1.dat", "cube_v2.dat", "cube_v3.dat"]

    try:
        convert_velocity_field(args.input_file, output_files)
        print("\nConversion completed successfully!")
        
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main() 
