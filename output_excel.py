#########################################################
####      Functions to Create Excel Spreadsheets     ####
#########################################################
import pandas as pd # also required openpyxl to be installed for verbose features

def excel_final(df_input, df1, df2, df3, df4, df5, df6, testname, save_dir):
    """Output true and engineering final data to Excel spreadsheet."""
    with pd.ExcelWriter(save_dir + testname + '_T.xlsx') as book_T:
        df_input.to_excel(book_T, sheet_name='input_parameters', index=False)#, header=False)
        df1.to_excel(book_T, sheet_name='1-wave', index=False)
        df2.to_excel(book_T, sheet_name='2-wave', index=False)
        df3.to_excel(book_T, sheet_name='3-wave', index=False)
    #
    with pd.ExcelWriter(save_dir + testname + '_E.xlsx') as book_E:
        df_input.to_excel(book_E, sheet_name='input_parameters', index=False, header=False)
        df4.to_excel(book_E, sheet_name='1-wave', index=False)
        df5.to_excel(book_E, sheet_name='2-wave', index=False)
        df6.to_excel(book_E, sheet_name='3-wave', index=False)
    return

if __name__ == '__main__':
    print("output_excel.py requires import from HopPyBar.py")