import argparse
import pubchempy as pcp
import pandas as pd
import os

def getcid(compound_name):
    #remove duplicted
    rddf=compound_name.drop_duplicates()
    rddf=rddf.reset_index(drop=True)

    dict={'cid':[], 'iupac_name':[], 'molecular_weight':[],'molecular_formula':[],'tip':[]};cnt=0
    for i in rddf['name']:
        cnt+=1

        cs = pcp.get_compounds(i, 'name')
        if len(cs)==1:
            dict['cid'].append(cs[0].cid)
            dict['iupac_name'].append(cs[0].iupac_name)
            dict['molecular_weight'].append(cs[0].molecular_weight)
            dict['molecular_formula'].append(cs[0].molecular_formula)
            dict['tip'].append('allright')
        elif len(cs)>1:
            dict['cid'].append(cs[0].cid)
            dict['iupac_name'].append(cs[0].iupac_name)
            dict['molecular_weight'].append(cs[0].molecular_weight)
            dict['molecular_formula'].append(cs[0].molecular_formula)
            dict['tip'].append('multi_id')
        elif len(cs) < 1:
            dict['cid'].append('na')
            dict['iupac_name'].append('na')
            dict['molecular_weight'].append('na')
            dict['molecular_formula'].append('na')
            dict['tip'].append('na')

        if cnt%20==0:
            print(f'========{cnt} compounds finished, wait!==========')
            rddf1 = pd.concat([rddf, pd.DataFrame(dict)], axis=1)
            result = pd.merge(compound_name, rddf1)
            result.to_csv('compounds.csv', index=False)
            print('=========compounds.csv temporarily saved======')
    rddf1 = pd.concat([rddf, pd.DataFrame(dict)], axis=1)
    return(rddf1)

def getname(compound_name):
    #remove duplicted
    rddf=compound_name.drop_duplicates()
    rddf=rddf.reset_index(drop=True)

    dict={'iupac_name':[], 'molecular_weight':[],'molecular_formula':[]};cnt=0
    for i in rddf['cid']:
        cnt+=1
        if i:
            cs = pcp.Compound.from_cid(int(i))

            dict['iupac_name'].append(cs.iupac_name)
            dict['molecular_weight'].append(cs.molecular_weight)
            dict['molecular_formula'].append(cs.molecular_formula)
        else:
            dict['iupac_name'].append(None)
            dict['molecular_weight'].append(None)
            dict['molecular_formula'].append(None)
        if cnt%20==0:
            print(f'========{cnt} compounds finished, wait!==========')
            rddf1 = pd.concat([rddf, pd.DataFrame(dict)], axis=1)
            result = pd.merge(compound_name, rddf1)
            result.to_csv('compounds.csv', index=False)
            print('=========compounds.csv temporarily saved======')
    rddf1 = pd.concat([rddf, pd.DataFrame(dict)], axis=1)
    return(rddf1)

def read_input_file(input_file):
    if input_file.endswith('.xlsx'):
        data = pd.read_excel(input_file, header=None)
    elif input_file.endswith('.csv'):
        data = pd.read_csv(input_file, header=None)
    elif input_file.endswith('.txt'):
        data = pd.read_csv(input_file, header=None, delimiter='\t')
    else:
        raise ValueError("Unsupported file format")
    
    return data

def main():
    parser = argparse.ArgumentParser(description='获取化合物信息的程序')
    parser.add_argument('input', type=str, help='输入文件路径')
    parser.add_argument('-o', '--output', type=str, help='输出文件名')
    parser.add_argument('-m', '--mode', choices=['name', 'cid'], default='name', help='输入数据模式')

    args = parser.parse_args()

    input_file = args.input
    output = args.output if args.output else 'compounds.csv'
    if output and os.path.isdir(output):
        output_file = os.path.join(output, 'compounds.csv')
    elif output:
        output_file = output
    else:
        output_file = 'compounds.csv'
        
    mode = args.mode

    print('===========read the file=========')
    data = read_input_file(input_file)
    deduplicated_data = data.drop_duplicates()
    print('===========start getting infoamation=========')

    print(f'total {data.shape[0]} rows')
    print(f"After removing duplicates and nan, total {deduplicated_data.shape[0]} coumpounds")
    
    if mode == 'name':
        deduplicated_data.rename(columns={0: 'name'}, inplace=True)
        data.rename(columns={0: 'name'}, inplace=True)
        rd=getcid(deduplicated_data)
    elif mode == 'cid':
        deduplicated_data.rename(columns={0: 'cid'}, inplace=True)
        data.rename(columns={0: 'cid'}, inplace=True)
        rd=getname(deduplicated_data)
    
    #result = pd.concat([data, rd], axis=1)
    result = pd.merge(data, rd, how='left', on=[mode])
    result.to_csv(output_file, index=False)
    print(f'=========All done, see {output_file}\!=======')

if __name__ == '__main__':
    main()

#cs = pcp.get_compounds('Hydroxy-eicosatrienoic acid', 'name')
