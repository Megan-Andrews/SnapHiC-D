
import cooler



def coarsen(inputfile, outputfile, chrom, factor):
    # c = cooler.Cooler(inputfile)
    # m = c.matrix(balance=False, as_pixels=True).fetch(chrom)
    base_uri = inputfile
    output_uri = outputfile
    chunksize = 10000
    cooler.coarsen_cooler(base_uri, output_uri, factor, chunksize)

def main():
    chrom = 'chr22'
    resolution = 100000
    curr_resolution = 10000

    factor = resolution/curr_resolution

    inputdir = "/project/compbio-lab/scHi-C/Lee2019/Human_single_cell_10kb_cool/"
    coarsedir = "/project/compbio-lab/scHi-C/Lee2019/Human_single_cell_100kb_cool/"
    outdir = "/project/compbio-lab/scHi-C/Lee2019/100kb_imputed_cool/ODC_100kb_imputed_cool/"

    file_list = "/home/maa160/SnapHiC-D/experiments/2023-06-30/SnapHiC-D/file_lists/ODC_file_listA.txt"
    
    file_no = 1
    with open(file_list, "r") as file:
    # Iterate through each line in the file
        for line in file:
            print(file_no, line.strip())
            run(filename = line.strip(), chrom=chrom, resolution=resolution, factor= factor, inputdir=inputdir, coarsedir=coarsedir, outdir=outdir)
            file_no += 1

        


def run(filename, chrom, resolution, factor, inputdir, coarsedir, outdir):
    
    inputfile = inputdir + filename
    coarsefile = coarsedir + filename.replace("_10kb_", "_100kb_")
    # outputfile = outdir + filename.replace("_10kb_", "_100kb_").replace(".cool","_imputed.cool")
   
    # if not os.path.exists(os.path.join(coarsedir, coarsefile)):
    coarsen(inputfile, coarsefile, chrom, factor)
    # else:
        # print("already coarsened")



if __name__ == "__main__":
    main() 