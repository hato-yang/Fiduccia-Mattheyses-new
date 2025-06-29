#include "../include/main.h"

// This is a collection of objects that can be stored in the doubly-linked list

// ############################################################

// Integer_data is mostly a test data object for debugging
void initialize_integer(struct integer *new_integer, int data)
{
	new_integer->data = data;
}

void print_integer(struct integer *object)
{
	printf("%d ", object->data);
}

// ############################################################

// The most fundamental aspect of FM is the movement of cells to the correct partition
void initialize_cell(struct cell *new_cell, int identifier, int area, int area2)
{
	new_cell->identifier = identifier;
	new_cell->area = area;
	new_cell->area2 = area2;
	struct dll *nets = (struct dll *)malloc(sizeof(struct dll));
	initialize_dll(nets);
	new_cell->nets = nets;

	new_cell->partition = 0;
	new_cell->gain = 0;

	new_cell->cell_state = FREE;
}

// To find number of pins, call size_dll(cell->nets)
int dll_size(struct dll *list)
{
	return list->size;
}
int cell_pinsize(struct cell *cell)
{
	// Should be an O(1) call
	dll_size(cell->nets);
}

// Adds a net to the the beginning of the cell's netlist
void add_net_to_cell(struct cell *cell, struct net *net)
{
	insert_node(cell->nets, 0, net);
}

void print_cell(struct cell *cell)
{
	printf("%d", cell->identifier + 1);
}

// Consider deleting
void print_cell_area(struct cell *object)
{
	printf("%d ", object->area);
}

// Garbage collection for cells!
// Need to demalloc nets dll, then cell struct

// This is a many-to-many
void delete_cell(struct cell *cell)
{
	garbage_collection_dll(cell->nets, DO_NOT_DEALLOC_DATA);
	free(cell);
}

// ############################################################

void initialize_net(struct net *net, int identifier)
{
	net->identifier = identifier;

	net->free_cells = malloc(sizeof(struct dll));
	initialize_dll(net->free_cells);

	net->locked_cells = malloc(sizeof(struct dll));
	initialize_dll(net->locked_cells);

	int *num_cells_in_ = malloc(2 * sizeof(int));
	num_cells_in_[PARTITION_A] = 0;
	num_cells_in_[PARTITION_B] = 0;
	net->num_cells_in_ = num_cells_in_;

	net->number_of_cells = 0;
}

void print_net(struct net *net)
{
	struct dll *cells = net->free_cells;
	printf("Net %d has %d cell[s]\n", net->identifier, cells->size);
	print_dll(cells, CELL);
}

// Deletes the net and the metadata, but doesn't delete the cells
// Does not remove the references in NET_dll, NET_array -- that should be done manually
//   Avoids freeing memory twice
void delete_net(struct net *undesired_net)
{
	// Nets have a list of free cells and a list of locked cells, both of which need to be dealt with

	// Delete list of free cells
	struct dll *free_cells = undesired_net->free_cells;
	struct node *temp_node_for_cell = ((struct node *)free_cells->head)->next;
	delete_net_helper(undesired_net, free_cells, temp_node_for_cell);
	garbage_collection_dll(undesired_net->free_cells, DO_NOT_DEALLOC_DATA);

	struct dll *locked_cells = undesired_net->locked_cells;
	temp_node_for_cell = ((struct node *)locked_cells->head)->next;
	delete_net_helper(undesired_net, locked_cells, temp_node_for_cell);
	garbage_collection_dll(undesired_net->locked_cells, DO_NOT_DEALLOC_DATA);

	free(undesired_net->num_cells_in_);
	free(undesired_net);
}

void delete_net_helper(struct net *undesired_net, struct dll *cellist, struct node *temp_node_for_cell)
{
	// The cell data structure being accessed at each node
	struct cell *temp_cell;
	// The cell's netlist being accessed at each node
	struct dll *temp_cell_netlist;
	// Access the netlist
	struct node *temp_node_for_net;
	// temp_net for checking whether the net is in the netlist to be deleted
	struct net *temp_net;

	// Remove references to the net from every cell's netlist
	while (temp_node_for_cell != cellist->tail)
	{
		temp_cell = temp_node_for_cell->data_structure;
		temp_cell_netlist = temp_cell->nets;
		temp_node_for_net = ((struct node *)temp_cell_netlist->head)->next;
		// print_dll(temp_cell_netlist, NET);
		// Look through each net in the netlist
		while (temp_node_for_net != temp_cell_netlist->tail)
		{
			temp_net = temp_node_for_net->data_structure;
			if (temp_net == undesired_net)
			{
				// Take the node out of the node_list
				remove_node(temp_node_for_net, temp_cell_netlist);
				// Free the node (but not the net)
				free(temp_node_for_net);
				// We've removed the net from the cell's netlist, we can skip the rest of the netlist and move to the next cell that has the undesired_net
				break;
			}
			temp_node_for_net = temp_node_for_net->next;
		}
		// Move to next
		temp_node_for_cell = temp_node_for_cell->next;
		// print_dll(temp_cell_netlist, NET);
	}
}

// ############################################################

// Two partitions and a metadata struct get malloc'ed
void initialize_two_partitions(struct condensed *information)
{

	// This value is used to malloc the gain table in partitions
	int max_nets = calculate_max_nets_on_cell(information->CELL_array, information->CELL_array_size);
	information->max_nets = max_nets;

	struct partition *partition_A = malloc(sizeof(struct partition));
	initialize_partition(partition_A, max_nets, PARTITION_A);
	struct partition *partition_B = malloc(sizeof(struct partition));
	initialize_partition(partition_B, max_nets, PARTITION_B);

	struct partition **access_ = malloc(2 * sizeof(struct partition *));
	access_[PARTITION_A] = partition_A;
	access_[PARTITION_B] = partition_B;
	information->access_ = access_;

	information->partition_A = partition_A;
	information->partition_B = partition_B;
}

// Useful for calculating the gain table size
int calculate_max_nets_on_cell(struct cell **CELL_array, int CELL_array_size)
{
	int max_nets = 0;
	int i, num_nets;
	for (i = 0; i < CELL_array_size; i++)
	{
		num_nets = CELL_array[i]->nets->size;
		if (num_nets > max_nets)
			max_nets = num_nets;
	}
	return max_nets;
}

void initialize_partition(struct partition *partition, int max_nets, int which_partition)
{
	// Create main gain array, Gain from [-n, n)
	int k = 2;																  ///////gain擴大2倍
	struct dll **GAIN_array = calloc(2 * k * max_nets, sizeof(struct dll *)); // 2 * max_nets, 修改成 2 * k *  max_nets,
	partition->GAIN_array = GAIN_array;
	partition->GAIN_array_size = 2 * k * max_nets; // 2 * max_nets, 修改成 2 * k *  max_nets,
	// printf("max_nets(initialize_partition):%d\n", max_nets);
	// printf("GAIN_array_size:%d\n", partition->GAIN_array_size);
	int i;
	struct dll *list_of_cells_with_same_gain;
	for (i = 0; i < 2 * k * max_nets; i++) // 2 * max_nets, 修改成 2 * k *  max_nets
	{
		list_of_cells_with_same_gain = malloc(sizeof(struct dll));
		initialize_dll(list_of_cells_with_same_gain);
		GAIN_array[i] = list_of_cells_with_same_gain;
	}

	partition->which_partition = which_partition;

	// Create list of cells
	struct dll *cells_in_partition = malloc(sizeof(struct dll));
	initialize_dll(cells_in_partition);
	partition->cells_in_partition = cells_in_partition;
	// Partition starts with 0 cell area
	partition->total_partition_area = 0;
	// Set the pointer to null
	partition->max_gain_pointer = NULL;
}

void populate_partitions(struct condensed *information)
{
	// This should be determined by the PARTITION_CELLS option in main.h
	if (PARTITION_CELLS == GENETIC_ALGORITHM)
		segregate_cells_with_GA(information);
	else
		segregate_cells_randomly(information); // segregate_cells_by_net_order(information); //

	int current_cutstate = calculate_initial_cutstate(information->NET_array, information->NET_array_size, information);

	// Generate cutstate list
	if (information->FM_chromosome == NULL)
		information->lowest_cutstate = current_cutstate;
}

void delete_partition(struct partition *undesired_partition)
{
	// Straight list of cells
	garbage_collection_dll(undesired_partition->cells_in_partition, DO_NOT_DEALLOC_DATA);
	int i;
	for (i = 0; i < undesired_partition->GAIN_array_size; i++)
	{
		garbage_collection_dll(undesired_partition->GAIN_array[i], DO_NOT_DEALLOC_DATA);
	}
	free(undesired_partition->GAIN_array);
	free(undesired_partition);
}
void populate_partitions_from_chromosome(struct condensed *information) // FM_REPEAT使用
{
	int i;
	for (i = 0; i < information->CELL_array_size; i++)
	{
		struct cell *c = information->CELL_array[i];
		int target_partition = information->FM_chromosome->gene_array[i];

		c->which_partition = target_partition;
		c->partition = information->access_[target_partition];

		if (target_partition == PARTITION_A)
		{
			insert_node(information->partition_A->cells_in_partition, 0, c);
			information->partition_A->total_partition_area += c->area;
		}
		else
		{
			insert_node(information->partition_B->cells_in_partition, 0, c);
			information->partition_B->total_partition_area += c->area2;
		}

		// 更新 net 中 cell 的分區計數
		struct node *net_node = c->nets->head->next;
		while (net_node != c->nets->tail)
		{
			struct net *net = (struct net *)net_node->data_structure;
			net->num_cells_in_[target_partition]++;
			insert_node(net->free_cells, 0, c);
			net_node = net_node->next;
		}
	}
}