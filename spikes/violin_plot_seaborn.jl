using PyCall, PyPlot
# Use Conda.jl to install seaborn unless already installed: using Conda; Conda.add("seaborn")
@pyimport seaborn as sns
sns.set_style("whitegrid")
tips = sns.load_dataset("tips")
sns.violinplot(x="day", y="total_bill", hue="smoker", data=tips, palette="muted", split=true)
filename = "violin_plot_example.png"
savefig(filename)
println("Saved to file: ", filename)