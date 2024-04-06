"""
Utils to write input.txt and coord_file for TEM-Simulator.
"""
def write_inputfile(base_address = "/content/cryoEM-notebooks/notebooks/Material/Data/TMV/", 
                    log_file = "simulator.log",
                    sample_diameter = 1000,
                    thickness_edge = 100,
                    thickness_center = 50,
                    pdb_file_in = "20M3.pdb",
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
  f.write("use_quantization = no\n")
  f.write("dqe = 0.8\n")
  f.write("mtf_a = 0.7\n")
  f.write("mtf_b = 0.2\n")
  f.write("mtf_c = 0.1\n")
  f.write("mtf_alpha = 10\n")
  f.write("mtf_beta = 40\n")
  f.write(f"image_file_out = {base_address+image_file_out}\n")
  f.write("=== detector ===\n")
  f.write("det_pix_x = 400\n")
  f.write("det_pix_y = 400\n")
  f.write("pixel_size = 15\n")
  f.write("gain = 10\n")
  f.write("use_quantization = no\n")
  f.write("dqe = 0.8\n")
  f.write("mtf_a = 0.7\n")
  f.write("mtf_b = 0.2\n")
  f.write("mtf_c = 0.1\n")
  f.write("mtf_alpha = 10\n")
  f.write("mtf_beta = 40\n")
  f.write(f"image_file_out = {base_address+image_file_out_nonoise}")
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
