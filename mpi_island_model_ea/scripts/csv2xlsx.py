#!/usr/bin/env python

import os
import pandas as pd
import shutil
import re
import openpyxl

head = os.getcwd()
extensions = ('.csv')

pd.set_option('display.max_colwidth', 90)

for subdir, dirs, files in os.walk(head):
    #tree = re.match(r'\/?([^\/]+)', relp).groups()
    #tree = re.match(r'\/?([^\/]+)', relp)
    for file in files:
        ext = os.path.splitext(file)[-1].lower()
        if ext in extensions and not file.startswith('.') and 'run' not in file:
            relp = os.path.relpath(subdir, head)
            tree = relp.split("/")
            combo = 'all_' + tree[0] + '.xlsx';
            if os.path.exists(combo):
                writer = pd.ExcelWriter(combo, engine='openpyxl', mode='a')
            else: 
                writer = pd.ExcelWriter(combo, engine='openpyxl')
                s1 = writer.book.create_sheet(title='sol')
                s2 = writer.book.create_sheet(title='top')
            print('found matching file: ' + file)
            fn = os.path.join(os.path.relpath(subdir, head), file)
            fout = fn + '.xlsx'
            dcp = head + '/all'
            fcp = dcp + '/' + os.path.split(fout)[1]
            if not os.path.isdir(dcp):
                os.mkdir(dcp, mode = 0o740)
            print('reading: ' + fn)
            df = pd.read_csv(fn)
            #print(df)
            sh = tree[0] + '-' + tree[1]
            print('worksheet: ' + sh)
            df.to_excel(writer, sheet_name=sh, header=True, index=False)
            writer.sheets = dict((ws.title, ws) for ws in writer.book.worksheets)
            if s1.max_column == 1:
                sol1 = df.iloc[:,[0]]
                sol1.to_excel(writer, sheet_name='sol', header=True, index=False)
                sol1.to_excel(writer, sheet_name='top', header=True, index=False)
            sol1 = df.iloc[:,[2,4]]
            sol2 = df.iloc[:,[12,19]]
            if 'bench' in file:
                sol1.columns = [ 'bench_avg_sol_fitness', 'bench_global_best_sol_fitness' ]
                sol2.columns = [ 'bench_avg_topo_fitness', 'bench_global_best_topo_fitness' ]
            if 'evo' in file:
                sol1.columns = [ 'evo_avg_sol_fitness', 'evo_global_best_sol_fitness' ]
                sol2.columns = [ 'evo_avg_topo_fitness', 'evo_global_best_topo_fitness' ]
            sol1.to_excel(writer, sheet_name='sol', startcol=s1.max_column, header=True, index=False)
            sol2.to_excel(writer, sheet_name='top', startcol=s2.max_column, header=True, index=False)
            for name, ws in writer.sheets.items():
                for column in ws.columns:
                    ws.column_dimensions[column[0].column_letter].width = 20
            writer.save()


