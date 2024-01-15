import pandas as pd
import numpy as np


class Stats:

    def __init__(self, df):
        self.df = df

    def add_mean(self, data,df):
        for value in df.mean().values:
            data.append(value)
        return data

    def add_std(self, data, df):
        for value in df.std().values:
            data.append(value)
        return data

    def add_cov(self, data_cov, data_mean, data_std):
        for i in range(1, len(data_mean)):
            data_cov.append(100 * (data_std[i] / data_mean[i]))
        return data_cov

    def add_mean_std_cov(self, tangential_exps=False, coronal_exps=False, all_exps=True, output_filename=None):
        new_df = pd.DataFrame(columns=self.df.columns)
        if tangential_exps:
            # Add mean and std
            tangential_inds_1 = np.logical_or(np.array(self.df['Exp_Name'] == 'MG49_lhs'),\
                                             np.array(self.df['Exp_Name'] == 'MG49_rhs'),)
            tangential_inds_2 = np.logical_or(np.array(self.df['Exp_Name'] == 'MG50_lhs'),\
                                             np.array(self.df['Exp_Name'] == 'MG50_rhs'))

            tangential_df = self.df[np.logical_or(tangential_inds_1,tangential_inds_2)]
            data_mean = ['Tangential_Mean']
            data_mean = self.add_mean(data_mean,tangential_df)
            new_df = new_df.append(pd.DataFrame(columns=new_df.columns, data=[data_mean]))
            data_std = ['Tangential_STD']
            data_std = self.add_std(data_std,tangential_df)
            new_df = new_df.append(pd.DataFrame(columns=new_df.columns, data=[data_std]))
            data_cov = ['Tangential_Coefficient_Of_Variation_%']
            data_cov = self.add_cov(data_cov, data_mean, data_std)
            new_df = new_df.append(pd.DataFrame(columns=new_df.columns, data=[data_cov]))

        if coronal_exps:
            coronal_inds = np.logical_or(np.array(self.df['Exp_Name'] == 'MG48_lhs'),\
                                         np.array(self.df['Exp_Name'] == 'MG48_rhs'))
            coronal_df = self.df[coronal_inds]
            data_mean = ['Coronal_Mean']
            data_mean = self.add_mean(data_mean,coronal_df)
            new_df = new_df.append(pd.DataFrame(columns=new_df.columns, data=[data_mean]))
            data_std = ['Coronal_STD']
            data_std = self.add_std(data_std,coronal_df)
            new_df = new_df.append(pd.DataFrame(columns=new_df.columns, data=[data_std]))
            data_cov = ['Coronal_Coefficient_Of_Variation_%']
            data_cov = self.add_cov(data_cov, data_mean, data_std)
            new_df = new_df.append(pd.DataFrame(columns=new_df.columns, data=[data_cov]))

        if all_exps:
            data_mean = ['Mean']
            data_mean = self.add_mean(data_mean,self.df)
            new_df = new_df.append(pd.DataFrame(columns=new_df.columns, data=[data_mean]))
            data_std = ['STD']
            data_std = self.add_std(data_std,self.df)
            new_df = new_df.append(pd.DataFrame(columns=new_df.columns, data=[data_std]))
            data_cov = ['Coefficient_Of_Variation_%']
            data_cov = self.add_cov(data_cov, data_mean, data_std)
            new_df = new_df.append(pd.DataFrame(columns=new_df.columns, data=[data_cov]))

        self.df = self.df.append(new_df)

        if output_filename is not None:
            self.df.to_csv(output_filename)

        return self.df
