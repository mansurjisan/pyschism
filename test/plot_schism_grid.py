from pyschism import read_gr3_file, plot_gr3

if __name__ == "__main__":
    filename = 'hgrid.gr3'
    output_file = 'schism_grid_plot.png'

    print(f"Reading file: {filename}")
    nodes, elements = read_gr3_file(filename)
    if nodes is not None and elements is not None:
#        plot_gr3(nodes, elements, "output_custom_mesh_only.png", min_lon=-86, max_lon=-78, min_lat=23, max_lat=30, show_land=True)
      plot_gr3(nodes, elements, "output_plot.png", min_lon=-88, max_lon=-77, min_lat=20, max_lat=30, show_land=False)

        
    else:
        print("Failed to read the .gr3 file. Please check the file path and format.")
