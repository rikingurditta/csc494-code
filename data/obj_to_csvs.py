from os.path import basename, splitext


def obj_to_csvs(file):
    v = []
    e = []
    for line in map(lambda x: x.strip('\n'), filter(lambda x: len(x.split(' ')) >= 3, file.readlines())):
        entries = line.split(' ')
        if entries[0] == 'v':
            v.append(f'{entries[1]},{entries[3]}')
        elif entries[0] == 'l':
            e.append(f'{int(entries[1]) - 1},{int(entries[2]) - 1}')
    with open(f'{splitext(basename(file.name))[0]}_V.csv', 'w') as v_file:
        v_file.write('\n'.join(v))
    with open(f'{splitext(basename(file.name))[0]}_E.csv', 'w') as e_file:
        e_file.write('\n'.join(e))


if __name__ == '__main__':
    files = ['2d_torus.obj',
             '2d_torus_dec1.obj',
             '2d_torus_dec2.obj',
             '2d_torus_dec3.obj']
    for filename in files:
        with open(filename) as f:
            obj_to_csvs(f)
