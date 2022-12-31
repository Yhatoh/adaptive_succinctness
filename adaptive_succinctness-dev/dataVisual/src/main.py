import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys


if __name__ == '__main__':
    input_file = sys.argv[1]
    df = pd.read_csv(input_file)
    output_file = input_file.replace('.csv', '')
    #print(output_file)
    df2 = df.sort_values(by=['total bpp'])
    #print(df)
    #print(df2)
    
    df2.to_excel(f'{output_file}_sorted.xlsx')
    df.to_excel(f'{output_file}.xlsx')