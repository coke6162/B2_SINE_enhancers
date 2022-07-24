#!/usr/bin/env python
# coding: utf-8

import plotly
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import os

# Read table
df = pd.read_csv('all_filtered.giggle.txt', sep='\t')
#print(df)

# Filter table to keep positive scores only (i.e. giggle_score>0)
columns = ['giggle_score']
filter_= (df[columns] > 0).all(axis=1)
positive_df = df[filter_]
#print(positive_df)

# Make a bubble plot
fig = px.scatter(positive_df, y="repeat_family", x="sample", size="giggle_score",  color="odds_ratio", size_max=25, opacity=1,
                 color_continuous_scale=[(0.0, '#FD9367'), (0.33, '#C3305D'), (0.67, '#782D65'), (1, '#432967')],
                 labels={"sample":"Sample", "giggle_score":"Enrichment", "repeat":"Repeats"},
                                  category_orders={"sample":["picc_WT_UT_H3K27ac", "picc_WT_IFNG_2h_H3K27ac", "picc_WT_UT_STAT1", "picc_WT_IFNG_2h_STAT1", "plat_WT_UT_STAT1", "plat_WT_IFNG_1.5h_STAT1"],
                                 "repeat_family":["B2_Mm2", "RLTR30B_MM", "RLTR30E_MM", "RLTR30D_MM", "PB1D7", "ID_B1", "L1MB2", "MT2B2", "MIR1_Amn", "MIR3", "LTR104_Mam"]}
                )

# Update plot layout
for template in ["plotly_white"]:    
    fig.update_layout(template=template)
    fig.update_layout(autosize=False,width=750,height=900)
    
fig.update_layout(font=dict(size=20))
fig.update_yaxes(title_standoff=20)
fig.update_xaxes(title_standoff=20)
fig.update_yaxes(tickfont_size=20, ticks="outside", tickangle=0, ticklen=4)
fig.update_xaxes(tickfont_size=20, tickangle=45, ticks="outside", ticklen=4)
fig.layout.font.family = 'Helvetica'

# Show and save fig
fig.show()
fig.write_image("giggle_bubbles.svg")