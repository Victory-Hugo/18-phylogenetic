# import pandas as pd
# import numpy as np

# # 读取数据
# file_path = 'C:/Users/victo/Desktop/Ancestor_Hap_Results.txt'
# output_path = 'C:/Users/victo/Desktop/Rho统计共祖时间.txt'
# df = pd.read_csv(file_path, sep='\t')

# # 需要的列名
# required_columns = [
#     'RHO_CS', 'RHO_CS_SE',
#     'RHO_SYN', 'RHO_SYN_SE',
#     'RHO_HVSI', 'RHO_HVSI_SE',
#     'RHO_HVSI_TRANSI', 'RHO_HVSI_TRANSI_SE',
#     'RHO_HVSII', 'RHO_HVSII_SE',
#     'RHO_CR', 'RHO_CR_SE'
# ]

# # 检查并添加缺失的列，值全部设为0
# for col in required_columns:
#     if col not in df.columns:
#         df[col] = 0

# # 定义函数
# def calculate_ages(df, rho_col, se_col, multiplier):
#     df[f'{rho_col}_AGE'] = df[rho_col] * multiplier
#     df[f'{rho_col}_AGE_95LB'] = (df[rho_col] - 1.96 * df[se_col]) * multiplier
#     df[f'{rho_col}_AGE_95HB'] = (df[rho_col] + 1.96 * df[se_col]) * multiplier

# # 定义计算
# cs_multiplier = 3624
# cs_base = -40.2789
# cs_exp = -0.0263
# df['CS_AGE'] = cs_multiplier * (((np.exp(-(np.exp((df['RHO_CS'] - cs_base) * cs_exp)))) * 0.4794) / 0.4794 * df['RHO_CS'])
# df['CS_AGE_95LB'] = cs_multiplier * (((np.exp(-(np.exp(((df['RHO_CS'] - 1.96 * df['RHO_CS_SE']) - cs_base) * cs_exp)))) * 0.4794) / 0.4794 * (df['RHO_CS'] - 1.96 * df['RHO_CS_SE']))
# df['CS_AGE_95HB'] = cs_multiplier * (((np.exp(-(np.exp(((df['RHO_CS'] + 1.96 * df['RHO_CS_SE']) - cs_base) * cs_exp)))) * 0.4794) / 0.4794 * (df['RHO_CS'] + 1.96 * df['RHO_CS_SE']))

# # 定义倍率
# multipliers = {
#     'RHO_SYN': 7872,
#     'RHO_HVSI': 16677,
#     'RHO_HVSI_TRANSI': 18845,
#     'RHO_HVSII': 22388,
#     'RHO_CR': 9058
# }

# # 循环计算 
# for rho_col, multiplier in multipliers.items():
#     calculate_ages(df, rho_col, f'{rho_col}_SE', multiplier)

# # 生成输出文件
# result_columns = [
#     'CS_AGE', 'CS_AGE_95LB', 'CS_AGE_95HB',
#     'RHO_SYN_AGE', 'RHO_SYN_AGE_95LB', 'RHO_SYN_AGE_95HB',
#     'RHO_HVSI_AGE', 'RHO_HVSI_AGE_95LB', 'RHO_HVSI_AGE_95HB',
#     'RHO_HVSI_TRANSI_AGE', 'RHO_HVSI_TRANSI_AGE_95LB', 'RHO_HVSI_TRANSI_AGE_95HB',
#     'RHO_HVSII_AGE', 'RHO_HVSII_AGE_95LB', 'RHO_HVSII_AGE_95HB',
#     'RHO_CR_AGE', 'RHO_CR_AGE_95LB', 'RHO_CR_AGE_95HB'
# ]
# result = df[result_columns]
# result.to_csv(output_path, sep='\t', encoding='UTF-8')

# # 打印
# print(result)
import pandas as pd
import numpy as np

# 读取数据
file_path = 'C:/Users/victo/Desktop/output_results.txt'
output_path = 'C:/Users/victo/Desktop/Rho统计共祖时间.txt'
df = pd.read_csv(file_path, sep='\t')

# 需要的列名
required_columns = [
    'Haplogroup', 'RHO_CS', 'RHO_CS_SE',
    'RHO_SYN', 'RHO_SYN_SE',
    'RHO_HVSI', 'RHO_HVSI_SE',
    'RHO_HVSI_TRANSI', 'RHO_HVSI_TRANSI_SE',
    'RHO_HVSII', 'RHO_HVSII_SE',
    'RHO_CR', 'RHO_CR_SE'
]

# 检查并添加缺失的列，值全部设为0
for col in required_columns:
    if col not in df.columns:
        df[col] = 0

# 定义函数
def calculate_ages(df, rho_col, se_col, multiplier):
    df[f'{rho_col}_AGE'] = df[rho_col] * multiplier
    df[f'{rho_col}_AGE_95LB'] = (df[rho_col] - 1.96 * df[se_col]) * multiplier
    df[f'{rho_col}_AGE_95HB'] = (df[rho_col] + 1.96 * df[se_col]) * multiplier

# 定义计算
cs_multiplier = 3624
cs_base = -40.2789
cs_exp = -0.0263
df['CS_AGE'] = cs_multiplier * (((np.exp(-(np.exp((df['RHO_CS'] - cs_base) * cs_exp)))) * 0.4794) / 0.4794 * df['RHO_CS'])
df['CS_AGE_95LB'] = cs_multiplier * (((np.exp(-(np.exp(((df['RHO_CS'] - 1.96 * df['RHO_CS_SE']) - cs_base) * cs_exp)))) * 0.4794) / 0.4794 * (df['RHO_CS'] - 1.96 * df['RHO_CS_SE']))
df['CS_AGE_95HB'] = cs_multiplier * (((np.exp(-(np.exp(((df['RHO_CS'] + 1.96 * df['RHO_CS_SE']) - cs_base) * cs_exp)))) * 0.4794) / 0.4794 * (df['RHO_CS'] + 1.96 * df['RHO_CS_SE']))

# 定义倍率
multipliers = {
    'RHO_SYN': 7872,
    'RHO_HVSI': 16677,
    'RHO_HVSI_TRANSI': 18845,
    'RHO_HVSII': 22388,
    'RHO_CR': 9058
}

# 循环计算 
for rho_col, multiplier in multipliers.items():
    calculate_ages(df, rho_col, f'{rho_col}_SE', multiplier)

# 生成输出文件
result_columns = [
    'Haplogroup', 'CS_AGE', 'CS_AGE_95LB', 'CS_AGE_95HB',
    'RHO_SYN_AGE', 'RHO_SYN_AGE_95LB', 'RHO_SYN_AGE_95HB',
    'RHO_HVSI_AGE', 'RHO_HVSI_AGE_95LB', 'RHO_HVSI_AGE_95HB',
    'RHO_HVSI_TRANSI_AGE', 'RHO_HVSI_TRANSI_AGE_95LB', 'RHO_HVSI_TRANSI_AGE_95HB',
    'RHO_HVSII_AGE', 'RHO_HVSII_AGE_95LB', 'RHO_HVSII_AGE_95HB',
    'RHO_CR_AGE', 'RHO_CR_AGE_95LB', 'RHO_CR_AGE_95HB'
]
result = df[result_columns]
result.to_csv(output_path, sep='\t', index=False, encoding='UTF-8')

# 打印
print(result)
