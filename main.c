#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct Node Node;

enum SpinState
{
    POS = 1,
    NEG = -1
};

struct Node
{
    int N_neighbours;
    Node *neighbours;
    int spin;
};

void display_matrix(Node **matrix, int side_length)
{
    for (size_t i = 0; i < side_length; i++)
    {
        for (size_t j = 0; j < side_length; j++)
        {
            if (matrix[i][j].spin == POS)
                printf(" %d  ", matrix[i][j].spin);
            else
                printf("%d  ", matrix[i][j].spin);
        }
        printf("\n");
    }
    printf("\n");
}

int metropolis_step(Node **matrix, int side_length, double J, double B)
{

    int rand_row = rint(rand() / (double)(RAND_MAX - 1) * (side_length - 1));
    int rand_col = rint(rand() / (double)(RAND_MAX - 1) * (side_length - 1));

    matrix[rand_col][rand_row].spin *= (-1); // Flip polarization of 1 random spin

    int c_spin = matrix[rand_col][rand_row].spin;
    int t_spin = matrix[(rand_col + 1) % side_length][rand_row].spin;
    int r_spin = matrix[rand_col][(rand_row + 1) % side_length].spin;
    int b_spin = matrix[(rand_col - 1 + side_length) % side_length][rand_row].spin;
    int l_spin = matrix[rand_col][(rand_row - 1 + side_length) % side_length].spin;

    double delta_energy = 0;
    delta_energy += 2 * J * (c_spin * t_spin > 0 ? 1 : -1);
    delta_energy += 2 * J * (c_spin * r_spin > 0 ? 1 : -1);
    delta_energy += 2 * J * (c_spin * b_spin > 0 ? 1 : -1);
    delta_energy += 2 * J * (c_spin * l_spin > 0 ? 1 : -1);

    delta_energy += -2 * B * c_spin;

    if (delta_energy < 0) return 0; // step accepted

    double a = rand() / (double)(RAND_MAX - 1);

    if (a > exp(delta_energy))
        return 0; // step accepted
    else
        matrix[rand_col][rand_row].spin *= (-1);
    return 0;
}

double compute_magnetization(Node **matrix, int side_length)
{
    double tmp_mag = 0;

    for (size_t i = 0; i < side_length; i++)
    {
        for (size_t j = 0; j < side_length; j++)
        {
            tmp_mag += matrix[i][j].spin;
        }
    }

    return tmp_mag / (side_length * side_length);
}

double compute_matrix_energy(Node **matrix, int side_length, double J, double B)
{

    double energy = 0;
    for (size_t i = 0; i < side_length; i++)
    {
        for (size_t j = 0; j < side_length; j++)
        {
            // Avoid double counting by sum only over the top and right neighbours
            energy += J * matrix[i][j].spin * matrix[(i + 1) % side_length][j].spin;
            energy += J * matrix[i][j].spin * matrix[i][(j + 1) % side_length].spin;

            energy += (-1) * B * matrix[i][j].spin;
        }
    }

    return energy;
}

int main(int argc, char const *argv[])
{

    srand(2); // Fix rand seed for reproducibility

    const int N_NODES_PER_SIDE = 10;
    double COUPLING_CONSTANT = 1;
    double EXTERNAL_MAGNETIC_FIELD = 3;

    // Node node_matrix[N_NODES_PER_SIDE][N_NODES_PER_SIDE];
    Node **node_matrix = (Node **)malloc(sizeof(Node *) * N_NODES_PER_SIDE);

    for (size_t i = 0; i < N_NODES_PER_SIDE; i++)
    {
        node_matrix[i] = (Node *)malloc(sizeof(Node) * N_NODES_PER_SIDE);
    }

    for (size_t i = 0; i < N_NODES_PER_SIDE; i++)
    {
        for (size_t j = 0; j < N_NODES_PER_SIDE; j++)
        {
            node_matrix[i][j] = (Node){
                .N_neighbours = 4,
                // 0 Top, 1 Right, 2 Bottom, 3 Left
                .neighbours = (Node *)malloc(sizeof(Node) * 4),
                // Random value -1 or 1
                .spin = rand() > (RAND_MAX / 2) ? POS : NEG,
            };

            node_matrix[i][j].neighbours[0] = node_matrix[(i + 1) % N_NODES_PER_SIDE][j];
            node_matrix[i][j].neighbours[1] = node_matrix[i][(j + 1) % N_NODES_PER_SIDE];
            node_matrix[i][j].neighbours[2] = node_matrix[(i - 1) % N_NODES_PER_SIDE][j];
            node_matrix[i][j].neighbours[3] = node_matrix[i][(j - 1) % N_NODES_PER_SIDE];
        }
    }

    display_matrix(node_matrix, N_NODES_PER_SIDE);

    for (size_t i = 0; i < 500; i++)
    {
        metropolis_step(node_matrix, N_NODES_PER_SIDE, COUPLING_CONSTANT, EXTERNAL_MAGNETIC_FIELD);
    }
    printf("%lf\n", compute_magnetization(node_matrix, N_NODES_PER_SIDE));

    printf("EXTERNAL REMOVED");
    EXTERNAL_MAGNETIC_FIELD = 0;

    for (size_t i = 0; i < 500; i++)
    {
        metropolis_step(node_matrix, N_NODES_PER_SIDE, COUPLING_CONSTANT, EXTERNAL_MAGNETIC_FIELD);
        printf("%lf\n", compute_magnetization(node_matrix, N_NODES_PER_SIDE));
    }

    display_matrix(node_matrix, N_NODES_PER_SIDE);

    return 0;
}
