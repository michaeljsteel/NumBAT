# User NumBAT mesh implementation file

nmtoum = 0.001  # template radii are in nm but matplotlib plots are in microns


#LINEWIDTH=0.75

#TODO:
# expand allowed_/reqd_ parameters to take a 'is_dim' field.
# have structure.Structure determine whether units are microns or nm.
# multiply by nm_to_um  in get_param if required.



def _process_one_and_two_incls(params):
    nelts = 0
    gmshfile = ""

    if "slab_b_x" in params:
        raise ValueError(
            f"NumBAT doesn't understand your geometry: with shape {params['inc_shape']}, I did not expect values for slab_b."
        )

    if "slab_a_x" in params:
        raise ValueError(
            f"NumBAT doesn't understand your geometry: with shape {params['inc_shape']}, I did not expect values for slab_a."
        )

    if "inc_a_x" in params:
        if (
            "coat_y" not in params and "inc_b_x" not in params
        ):  # One inclusion, no coating
            gmshfile = "oneincl"  # used to be just '1'
            nelts = 2  # bkg, core (mat_a)

        elif (
            "coat_y" not in params and "inc_b_x" in params
        ):  # Two inclusions, no coating
            gmshfile = "twoincl"  # used to be just '2'
            nelts = 3  # bkg, core 1 (mat_a), core 2 (mat_b)

        # Two inclusions, with coating # TODO:implement
        elif "coat_y" in params and "inc_b_x" in params:
            raise NotImplementedError("Have not implemented 2 coated inclusions.")

        elif (
            "coat_y" in params and "inc_b_x" not in params
        ):  # One inclusion, with coating # TODO:implement
            raise NotImplementedError("Have not implemented 1 coated inclusions.")

        else:
            raise ValueError("NumBAT doesn't understand your geometry.")
    else:
        raise ValueError("must have at least one nonzero inclusion.")

    return gmshfile, nelts


def _process_one_and_two_incls_subs(msh_template):
    # TODO: these are crazy small defaults
    subs = [
        ("dx_in_nm = 100;", "dx_in_nm = %f;", "domain_x"),
        ("dy_in_nm = 50;", "dy_in_nm = %f;", "domain_y"),
        ("a1 = 20;", "a1 = %f;", "inc_a_x"),
        ("a1y = 10;", "a1y = %f;", "inc_a_y"),
        ("lc = 0.1;", "lc = %f;", "lc_bkg"),
        ("lc_refine_1 = lc/1;", "lc_refine_1 = lc/%f;", "lc_refine_1"),
        ("lc_refine_2 = lc/1;", "lc_refine_2 = lc/%f;", "lc_refine_2"),
    ]

    if msh_template in ["twoincl", "2", "2_on_s", "2_on_2s"]:
        subs.append(("a2 = 10;", "a2 = %f;", "inc_b_x"))
        subs.append(("a2y = 20;", "a2y = %f;", "inc_b_y"))
        subs.append(("sep = 10;", "sep = %f;", "two_inc_sep"))

        # geo = geo.replace('lc_refine_3 = lc/1;', 'lc_refine_3 = lc/%f;' % umb.lc_refine_3)
    if msh_template == "2":
        subs.append(("yoff = -5;", "yoff = %f;", "incs_y_offset"))

    if msh_template in ["1_on_slab", "1_on_2slabs", "1_on_slab", "2_on_2slabs"]:
        subs.append(("slab_width = dx_in_nm;", "slab_width = %f;", "slab_a_x"))
        subs.append(("slab_height = 10;", "slab_height = %f;", ".slab_a_y"))
        subs.append(("lc_refine_3 = lc/1;", "lc_refine_3 = lc/%f;", "lc_refine_3"))
        subs.append(("lc_refine_4 = lc/1;", "lc_refine_4 = lc/%f;", "lc_refine_4"))

    if msh_template in ["1_on_2slabs", "2_on_2slabs"]:
        subs.append(("slab2_width = dx_in_nm;", "slab2_width = %f;", "slab_b_x"))
        subs.append(("slab2_height = 5;", "slab2_height = %f;", "slab_b_y"))
        subs.append(("lc_refine_3 = lc/1;", "lc_refine_3 = lc/%f;", ".lc_refine_3"))
        subs.append(("lc_refine_4 = lc/1;", "lc_refine_4 = lc/%f;", "lc_refine_4"))

    return subs
