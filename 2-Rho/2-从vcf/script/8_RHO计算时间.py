#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np

def calculate_ages(df, rho_col, se_col, multiplier):
    df[f'{rho_col}_AGE']     = df[rho_col] * multiplier
    df[f'{rho_col}_AGE_95LB'] = (df[rho_col] - 1.96 * df[se_col]) * multiplier
    df[f'{rho_col}_AGE_95HB'] = (df[rho_col] + 1.96 * df[se_col]) * multiplier

def main():
    parser = argparse.ArgumentParser(description="计算各 haplogroup 的共祖时间")
    parser.add_argument('-i', '--input',  required=True, help="输入 TSV 文件路径")
    parser.add_argument('-o', '--output', required=True, help="输出结果 TXT 路径")
    args = parser.parse_args()

    # 读取数据
    df = pd.read_csv(args.input, sep='\t')

    # 重命名第一列为 Haplogroup
    df.rename(columns={df.columns[0]: 'Haplogroup'}, inplace=True)

    # 确保所有必须的列都存在，不存在的补 0
    required_columns = [
        'Haplogroup', 'RHO_CS', 'RHO_CS_SE',
        'RHO_SYN', 'RHO_SYN_SE',
        'RHO_HVSI', 'RHO_HVSI_SE',
        'RHO_HVSI_TRANSI', 'RHO_HVSI_TRANSI_SE',
        'RHO_HVSII', 'RHO_HVSII_SE',
        'RHO_CR', 'RHO_CR_SE'
    ]
    for col in required_columns:
        if col not in df.columns:
            df[col] = 0

    # special CS age
    cs_multiplier = 3624
    cs_base       = -40.2789
    cs_exp        = -0.0263
    df['CS_AGE']     = cs_multiplier * (np.exp(-np.exp((df['RHO_CS'] - cs_base) * cs_exp)) * df['RHO_CS'])
    df['CS_AGE_95LB'] = cs_multiplier * (np.exp(-np.exp(((df['RHO_CS'] - 1.96 * df['RHO_CS_SE']) - cs_base) * cs_exp)) * 
                                        (df['RHO_CS'] - 1.96 * df['RHO_CS_SE']))
    df['CS_AGE_95HB'] = cs_multiplier * (np.exp(-np.exp(((df['RHO_CS'] + 1.96 * df['RHO_CS_SE']) - cs_base) * cs_exp)) * 
                                        (df['RHO_CS'] + 1.96 * df['RHO_CS_SE']))

    # 其他 rho 的倍率
    multipliers = {
        'RHO_SYN':          7872,
        'RHO_HVSI':        16677,
        'RHO_HVSI_TRANSI': 18845,
        'RHO_HVSII':       22388,
        'RHO_CR':           9058
    }
    for rho_col, multiplier in multipliers.items():
        calculate_ages(df, rho_col, f'{rho_col}_SE', multiplier)

    # 输出
    result_columns = [
        'Haplogroup', 'CS_AGE', 'CS_AGE_95LB', 'CS_AGE_95HB',
        'RHO_SYN_AGE', 'RHO_SYN_AGE_95LB', 'RHO_SYN_AGE_95HB',
        'RHO_HVSI_AGE', 'RHO_HVSI_AGE_95LB', 'RHO_HVSI_AGE_95HB',
        'RHO_HVSI_TRANSI_AGE', 'RHO_HVSI_TRANSI_AGE_95LB', 'RHO_HVSI_TRANSI_AGE_95HB',
        'RHO_HVSII_AGE', 'RHO_HVSII_AGE_95LB', 'RHO_HVSII_AGE_95HB',
        'RHO_CR_AGE', 'RHO_CR_AGE_95LB', 'RHO_CR_AGE_95HB'
    ]
    df[result_columns].to_csv(args.output, sep='\t', index=False, encoding='utf-8')
    print(df[result_columns])

if __name__ == '__main__':
    main()
