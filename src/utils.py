"""
Utils to write input.txt and coord_file for TEM-Simulator, run singleton TEM simulations, and multiplle sims in parallel.
"""

def run_parallel_TEM_simulations(args_array: np.array = None) -> np.array:
  """
  Takes a numpy array of dimesion [num_processes, 22]
  Runs num_processes TEM_simulations in parallel with the inputs.
  Returns a [num_processes, 400, 400] numpy array of images.
  """
  if args_array is None:
    args_array = np.array([[11, 200.0, 1.3, 6000, 30000.0, 2.0, 2.0, 50.0, 3.0, 0.1, 0.8, 0.7, 0.2, 0.1, 10.0, 40.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                       [12, 200.0, 1.3, 60000, 60000.0, 2.0, 2.0, 50.0, 3.0, 0.1, 0.8, 0.7, 0.2, 0.1, 10.0, 40.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                       ])
  sim_results = []
  with concurrent.futures.ProcessPoolExecutor() as executor:
    results = [executor.submit(run_TEM_simulation, args) for args in args_array]
  
    for f in concurrent.futures.as_completed(results):
      sim_results.append(f.result())

  return np.array(sim_results)



def run_TEM_simulation(args: np.array) -> np.array:
  """
  Takes a 22 dimensional vector of inputs: index, coords, input_file_args.
  Creates a folder for the data,
  Writes an input.txt and a coord.txt file in folder. Moves pdf files to the folder.
  Runs the TEM-Simulator on these files.
  Returns the micrograph of the simulation as an [1, dim_x, dim_y] numpy vector.
  Args:
    args: A Numpy vector of 22-dims,
    [index, acc_voltage, energy_spread, total_dose, magnification, cs, cc, 
    aperture, focal_length, cond_ap_angle, dqe, mtf_a, mtf_b, mtf_c, mtf_alpha, 
    mtf_beta, 
    x, y, z, theta_x, theta_y, theta_z]
  Returns:
    np.array of dimensions [1, 400, 400] (default dim_x = dim_y = 400)
  """
  # Unpack
  index, acc_voltage, energy_spread, total_dose, magnification, cs, cc, \
  aperture, focal_length, cond_ap_angle, dqe, mtf_a, mtf_b, mtf_c, mtf_alpha, \
  mtf_beta, x, y, z, theta_x, theta_y, theta_z = args.tolist()
  
  index = int(index)

  base_address = "/content/" + "simulation_" + str(index) + "/"
  if not os.path.exists(base_address):
    os.makedirs(base_address, exist_ok=True)
  #4 Files needed for sim: pdb, transfer, input and coords.
  #copy files for simulation to dir
    copy_tree("/content/TEM-Simulator", base_address+"TEM-Simulator")
    shutil.copyfile("/content/cryoEM-notebooks/notebooks/Material/Data/TMV/2OM3.pdb", base_address+"2OM3.pdb")
    shutil.copyfile("/content/cryoEM-notebooks/notebooks/Material/Data/TMV/2OM3_transf.txt", base_address+"2OM3_transf.txt")
  write_inputfiles(base_address=base_address,
                  acc_voltage=acc_voltage,
                  energy_spread = energy_spread,
                  total_dose=total_dose,
                  magnification=magnification,
                  cs=cs,
                  cc=cc,
                  aperture=aperture,
                  focal_length=focal_length,
                  cond_ap_angle=cond_ap_angle,
                  dqe=dqe,
                  mtf_a=mtf_a,
                  mtf_b=mtf_b,
                  mtf_c=mtf_c,
                  mtf_alpha=mtf_alpha,
                  mtf_beta=mtf_beta)
  write_coords(base_address, x_coord=x, y_coord=y, z_coord=z, x_theta=theta_x, y_theta=theta_y, z_theta=theta_z)
  os.system('{0} {1}'.format(base_address+"TEM-Simulator/TEM-Simulator", base_address+"input.txt"))
  data = cryoemio.mrc2data(mrc_file = base_address+"tiltseries.mrc")
  return data[0,...]

def write_inputfiles(base_address = "/content/cryoEM-notebooks/notebooks/Material/Data/TMV/", 
                    log_file = "simulator.log",
                    pdb_file_in = "2OM3.pdb",
                    pdb_transf_file_in = "2OM3_transf.txt",
                    map_file_re_out = "2OM3_map.mrc",
                    map_file_im_out = "2OM3_abs_map.mrc",
                    particle_coords = "file",
                    coord_file_in = "TMV_coord.txt",
                    acc_voltage = 200,
                    energy_spread = 1.3,
                    total_dose = 60000.0,
                    magnification = 30000.0,
                    cs = 2.0,
                    cc = 2.0,
                    aperture = 50.0,
                    focal_length = 3.0,
                    cond_ap_angle = 0.1,
                    dqe = 0.8,
                    mtf_a = 0.7,
                    mtf_b = 0.2,
                    mtf_c = 0.1,
                    mtf_alpha = 10.0,
                    mtf_beta = 40.0,
                    image_file_out="tiltseries.mrc",
                    ):
  fpath = base_address + "input.txt"
  f = open(fpath, "w")
  f.write("=== simulation ===\n")
  f.write("generate_micrographs = yes\n")
  f.write(f"log_file = {base_address+log_file}\n")
  f.write("=== sample ===\n")
  f.write(f"diameter = 1000\n")
  f.write(f"thickness_edge = 100\n")
  f.write(f"thickness_center = 50\n")
  f.write("=== particle TMV ===\n")
  f.write("source = pdb\n")
  f.write(f"pdb_file_in = {base_address+pdb_file_in}\n")
  f.write(f"pdb_transf_file_in = {base_address+pdb_transf_file_in}\n")
  f.write(f"voxel_size = 0.1\n")
  f.write(f"map_file_re_out = {base_address+map_file_re_out}\n")
  f.write(f"map_file_im_out = {base_address+map_file_im_out}\n")
  f.write("=== particleset ===\n")
  f.write(f"particle_type = TMV\n")
  f.write(f"num_particles = 1\n")
  f.write(f"particle_coords = {particle_coords}\n")
  f.write(f"coord_file_in = {base_address+coord_file_in}\n")
  f.write("=== geometry ===\n")
  f.write(f"gen_tilt_data = yes\n")
  f.write(f"ntilts = 1\n")
  f.write(f"theta_start = 0\n")
  f.write(f"theta_incr = 2\n")
  f.write(f"geom_errors = none\n")
  f.write("=== electronbeam ===\n")
  f.write(f"acc_voltage = {acc_voltage}\n")
  f.write(f"energy_spread = {energy_spread}\n")
  f.write(f"gen_dose = yes\n")
  f.write(f"total_dose = {total_dose}\n")
  f.write("=== optics ===\n")
  f.write(f"magnification = {magnification}\n")
  f.write(f"cs = {cs}\n")
  f.write(f"cc = {cc}\n")
  f.write(f"aperture = {aperture}\n")
  f.write(f"focal_length = {focal_length}\n")
  f.write(f"cond_ap_angle = {cond_ap_angle}\n")
  f.write("gen_defocus = yes\n")
  f.write("defocus_nominal = 5\n")
  f.write("=== detector ===\n")
  f.write("det_pix_x = 400\n")
  f.write("det_pix_y = 400\n")
  f.write("pixel_size = 15\n")
  f.write("gain = 10\n")
  f.write("use_quantization = yes\n")
  f.write(f"dqe = {dqe}\n")
  f.write(f"mtf_a = {mtf_a}\n")
  f.write(f"mtf_b = {mtf_b}\n")
  f.write(f"mtf_c = {mtf_c}\n")
  f.write(f"mtf_alpha = {mtf_alpha}\n")
  f.write(f"mtf_beta = {mtf_beta}\n")
  f.write(f"image_file_out = {base_address+image_file_out}\n")
  f.close()


def write_inputfile(base_address = "/content/cryoEM-notebooks/notebooks/Material/Data/TMV/", 
                    log_file = "simulator.log",
                    sample_diameter = 1000,
                    thickness_edge = 100,
                    thickness_center = 50,
                    pdb_file_in = "2OM3.pdb",
                    pdb_transf_file_in = "2OM3_transf.txt",
                    voxel_size = 0.1,
                    map_file_re_out = "2OM3_map.mrc",
                    map_file_im_out = "2OM3_abs_map.mrc",
                    particle_type = "TMV",
                    num_particles = 1,
                    particle_coords = "file",
                    coord_file_in = "TMV_coord.txt",
                    gen_tilt_data = "yes",
                    ntilts = 1,
                    theta_start = 0,
                    theta_incr = 2,
                    acc_voltage = 200,
                    energy_spread = 1.3,
                    image_file_out="tiltseries.mrc",
                    image_file_out_nonoise="tiltseries.mrc"
                    ):
  fpath = base_address + "input.txt"
  f = open(fpath, "w")
  f.write("=== simulation ===\n")
  f.write("generate_micrographs = yes\n")
  f.write(f"log_file = {base_address+log_file}\n")
  f.write("=== sample ===\n")
  f.write(f"diameter = {sample_diameter}\n")
  f.write(f"thickness_edge = {thickness_edge}\n")
  f.write(f"thickness_center = {thickness_center}\n")
  f.write("=== particle TMV ===\n")
  f.write("source = pdb\n")
  f.write(f"pdb_file_in = {base_address+pdb_file_in}\n")
  f.write(f"pdb_transf_file_in = {base_address+pdb_transf_file_in}\n")
  f.write(f"voxel_size = {voxel_size}\n")
  f.write(f"map_file_re_out = {base_address+map_file_re_out}\n")
  f.write(f"map_file_im_out = {base_address+map_file_im_out}\n")
  f.write("=== particleset ===\n")
  f.write(f"particle_type = {particle_type}\n")
  f.write(f"num_particles = {num_particles}\n")
  f.write(f"particle_coords = {particle_coords}\n")
  f.write(f"coord_file_in = {base_address+coord_file_in}\n")
  f.write("=== geometry ===\n")
  f.write(f"gen_tilt_data = {gen_tilt_data}\n")
  f.write(f"ntilts = {ntilts}\n")
  f.write(f"theta_start = {theta_start}\n")
  f.write(f"theta_incr = {theta_incr}\n")
  f.write(f"geom_errors = none\n")
  f.write("=== electronbeam ===\n")
  f.write(f"acc_voltage = {acc_voltage}\n")
  f.write(f"energy_spread = {energy_spread}\n")
  f.write(f"gen_dose = yes\n")
  f.write(f"total_dose = 6000\n")
  f.write("=== optics ===\n")
  f.write("magnification = 30000\n")
  f.write("cs = 2\n")
  f.write("cc = 2\n")
  f.write("aperture = 50\n")
  f.write("focal_length = 3\n")
  f.write("cond_ap_angle = 0.1\n")
  f.write("gen_defocus = yes\n")
  f.write("defocus_nominal = 5\n")
  f.write("=== detector ===\n")
  f.write("det_pix_x = 400\n")
  f.write("det_pix_y = 400\n")
  f.write("pixel_size = 15\n")
  f.write("gain = 10\n")
  f.write("use_quantization = yes\n")
  f.write("dqe = 0.8\n")
  f.write("mtf_a = 0.7\n")
  f.write("mtf_b = 0.2\n")
  f.write("mtf_c = 0.1\n")
  f.write("mtf_alpha = 10\n")
  f.write("mtf_beta = 40\n")
  f.write(f"image_file_out = {base_address+image_file_out}\n")
  f.close()


def write_coords(base_address: str = "/content/cryoEM-notebooks/notebooks/Material/Data/TMV/", 
                coord_file: str = "TMV_coord.txt",
                 num_particles: int = 1,
                 x_coord: int = 0,
                 y_coord: int = 0,
                 z_coord: int = 0,
                 x_theta: int = 45,
                 y_theta: int = 45,
                 z_theta: int = 45):
  """
  Writes coord_file for the position of the SINGLE particle to be simulated.
  Args:
    base_address: directory location to write file,
    coord_file: file to write,
    num_particles: Number of particles in simulation (Hard coded to 1),
    x_coord: x coordinate of particle center,
    y_coord: y coordinate of particle center,
    z_coord: z coordinate of particle center,
    x_theta: rotation of particle in degrees about x-axis,
    y_theta: rotation of particle in degrees about y-axis,
    z_theta: rotation of particle in degrees about z-axis.
  Returns:
    None
  """
  fpath = base_address + coord_file
  f = open(fpath, "w")
  f.write(f" {num_particles} 6\n")
  f.write(f"{x_coord} {y_coord} {z_coord} {x_theta} {y_theta} {z_theta}")
  f.close()
