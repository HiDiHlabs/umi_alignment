import os
import pandas as pd



column_set1=['Raw Reads Filtered Due to Insufficient Support',  'Raw Reads Filtered Due to Mismatching Cigars',  'Raw Reads Filtered Due to Low Base Quality',
    'Raw Reads Filtered Due to Orphan Consensus Created','Consensus reads emitted', 'Total Raw Reads Considered']
column_set2=['overlapping templates', 'corrected templates', 'overlapping bases',  'corrected bases' ]
consensus_df=pd.DataFrame(columns=column_set1+column_set2, index=[snakemake.params['sample']])
# print(consensus_df)
with open(snakemake.input[0]) as handle:
    for line in handle:
        for i in consensus_df.columns:
            if i in line:
                # print(line)
                if i in column_set1:
                    line=line.split(':')[-1].split('(')[0].replace('.', '').replace(',','').replace(' ', '').strip()
                else:
                    line=line.split(']')[-1].split('(')[0].replace('.', '').replace(i, '').replace(',','').replace(' ', '').strip()
                line=int(line)
                consensus_df.loc['T-NHL_98_malign', i]=line
                break

consensus_df.to_csv(snakemake.output[0], sep='\t', index=True)
