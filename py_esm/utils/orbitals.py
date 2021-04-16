ang_mom_characters = 'spdfghiklmnoqrtuvwxyzabce'

ang_mom_dict = {
    0: {
        0: (0, 0, 0)
    },
    1: {
        -1: (1, 0, 0),      # x
        0: (0, 1, 0),       # y
        1: (0, 0, 1),       # z
    },
    2: {
        -2: (1, 1, 0),      # xy
        -1: (1, 0, 1),      # xz
        0: (0, 0, 2),       # z^2
        1: (0, 1, 1),       # yz
        2: (0, 1, 1),       # x^2 - y^2 # todo fix this one
    },
    3: {
        -3: (),             # x(x^2 - 3y^2)
        -2: (1, 1, 1),      # xyz
        -1: (1, 0, 2),      # xz^2
        0: (0, 0, 3),       # z^3
        1: (0, 1, 2),       # yz^2
        2: (1, 1, 1),       # (x^2 - y^2)z # todo fix this one
        3: (),              # y(3x^2 - y^2)
    }
}