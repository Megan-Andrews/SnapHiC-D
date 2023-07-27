import pandas as pd
import os
import matplotlib.pyplot as plt

result_dir = "/project/compbio-lab/scHi-C/Lee2019/results/2023-07-25/"
batch_21_DCC = "/project/compbio-lab/scHi-C/Lee2019/simulations/Astro_MG_190315_21yr/DCC.txt"
batch_29_DCC = "/project/compbio-lab/scHi-C/Lee2019/simulations/Astro_MG_190315_29yr/DCC.txt"

batches = ["Astro_MG_190315_21yr", "Astro_MG_190315_29yr"]
experiment_sizes = [90, 45, 30, 18, 15, 10, 9, 6, 5, 3, 2]

filters = ["filterA", "filterB"]

batch_21_DCC = pd.read_csv(batch_21_DCC, sep="\t", header=None, names=["bin1_id", "bin2_id", "LogFC_ground_truth"])
batch_29_DCC = pd.read_csv(batch_29_DCC, sep="\t", header=None, names=["bin1_id", "bin2_id", "LogFC_ground_truth"])

def create_plots(batch, filter_type):
    print(batch, " ", filter_type)
    if batch == "Astro_MG_190315_21yr":
        batch_df = batch_21_DCC
    else:
        batch_df = batch_29_DCC

    result_df = pd.DataFrame(columns=["Error", "Experiment Size"])
    accuracy_df = pd.DataFrame(columns=["TP", "TN", "FP", "FN", "Experiment Size", "Accuracy"])

    for e in experiment_sizes:
        file_name = f"diffHiC_{batch}_{e}x{e}_{filter_type}_results.csv"
        file_path = os.path.join(result_dir, file_name)
        # Read the CSV file into a pandas DataFrame
        df = pd.read_csv(file_path)
        df = df[df["FDR"] < 0.05]
        print(df.columns)
        df = df[['start1', 'start2', 'logFC', 'FDR']]
        df = df.rename(columns={"start1": "bin1_id", "start2": "bin2_id"})
        df = df.merge(batch_df, on=["bin1_id", "bin2_id"], how="outer")
        df = df[['logFC', 'LogFC_ground_truth']]
        print(df)
        FN = df["logFC"].isna().sum()
        FP = df["LogFC_ground_truth"].isna().sum()
        TP = df.shape[0] - FN - FP
        TN = 514**2 - FP - FN - TP
        accuracy = (TP + TN) / (TP + TN + FP + FN)
        print(FN, FP, TP, TN, accuracy)
        # Create a new DataFrame with the row data
        new_row = pd.DataFrame({"TP": [TP], "TN": [TN], "FP": [FP], "FN": [FN], "Experiment Size": [e], "Accuracy": [accuracy]})

        # Concatenate the new DataFrame with the original DataFrame
        accuracy_df = pd.concat([accuracy_df, new_row], ignore_index=True)

        df.loc[df["LogFC_ground_truth"].notna() & df["logFC"].notna(), "accuracy"] = "TP"
        df.loc[df["LogFC_ground_truth"].isna() & df["logFC"].notna(), "accuracy"] = "FP"
        df.loc[df["LogFC_ground_truth"].notna() & df["logFC"].isna(), "accuracy"] = "FN"

        df = df.fillna(1)
        df["Error"] = (df["logFC"] - df["LogFC_ground_truth"])**2
        df["Experiment Size"] = e
        df = df[[ "accuracy","Error", "Experiment Size"]]

        result_df = pd.concat([result_df, df], ignore_index=True)

    # Plotting the data
    plt.figure()
    plt.scatter(result_df[result_df["accuracy"]=="TP"]["Experiment Size"], result_df[result_df["accuracy"]=="TP"]["Error"], alpha=0.7, c="blue")
    plt.xlabel("Experiment Size")
    plt.ylabel("LogFC Error^2")
    plt.title(f"LogFC Error - {batch} - {filter_type}")
    plt.savefig(f"TP_LFC_ground_{batch}_{filter_type}_plot.png")
    plt.close()

    plt.figure()
    plt.scatter(result_df[result_df["accuracy"]=="FP"]["Experiment Size"], result_df[result_df["accuracy"]=="FP"]["Error"], alpha=0.7, c="red")
    plt.xlabel("Experiment Size")
    plt.ylabel("LogFC Error^2")
    plt.title(f"LogFC Error - {batch} - {filter_type}")
    plt.savefig(f"FP_LFC_ground_{batch}_{filter_type}_plot.png")
    plt.close()

    plt.figure()
    plt.scatter(result_df[result_df["accuracy"]=="FN"]["Experiment Size"], result_df[result_df["accuracy"]=="FN"]["Error"], alpha=0.7, c="green")
    plt.xlabel("Experiment Size")
    plt.ylabel("LogFC Error^2")
    plt.title(f"LogFC Error - {batch} - {filter_type}")
    plt.savefig(f"FN_LFC_ground_{batch}_{filter_type}_plot.png")
    plt.close()

    plt.figure()
    plt.scatter(accuracy_df["Experiment Size"], accuracy_df["Accuracy"], alpha=0.7)
    plt.xlabel("Experiment Size")
    plt.ylabel("Accuracy")
    plt.title(f"Accuracy - {batch} - {filter_type}")
    plt.savefig(f"Accuracy_ground_{batch}_{filter_type}_plot.png")
    plt.close()


create_plots("Astro_MG_190315_21yr", "filterA")
create_plots("Astro_MG_190315_21yr", "filterB")
create_plots("Astro_MG_190315_29yr", "filterA")
create_plots("Astro_MG_190315_29yr", "filterB")

