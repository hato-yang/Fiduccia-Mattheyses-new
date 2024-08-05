#include "../include/main.h"

/////////////
/*
void print_partition_cell_ids(struct partition *partition)
{
	// FILE *da = fopen("partition", "a");
	printf("partition %d 's Cell ID: ", partition->which_partition);
	// fprintf(da, "partition %d 's Cell ID: \n", partition->which_partition);
	struct node *current_node = partition->cells_in_partition->head->next;
	while (current_node != partition->cells_in_partition->tail)
	{
		struct cell *cell = (struct cell *)current_node->data_structure;
		printf("%d ", cell->identifier); //
		// fprintf(da, "%d \n", cell->identifier);
		current_node = current_node->next;
	}
	// fclose(da);
	printf("\n");
}*/
//////////////
void fiduccia_mattheyses_algorithm(struct condensed *information)
{

	calculate_initial_gains_wrapper(information); // gain越高cut越多，優先移動，以?儲存?????

	int cells_can_still_be_moved = 1;
	int timestep = 0;
	///////////////
	/////////////////計算初始分區面積
	int initial_partition_areas_A = 0;
	int initial_partition_areas_B = 0;
	int area0;
	int area1;
	struct node *node0 = information->partition_A->cells_in_partition->head->next;
	while (node0 != information->partition_A->cells_in_partition->tail)
	{
		struct cell *cell0 = (struct cell *)node0->data_structure;
		area0 = 0;
		if (cell0->which_partition == 0)
			area0 = cell0->area;
		printf("%d      %d           %4d    %4d\n", cell0->identifier + 1, cell0->which_partition, cell0->area, cell0->area2); //
		initial_partition_areas_A = initial_partition_areas_A + area0;
		node0 = node0->next;
	}
	struct node *node1 = information->partition_B->cells_in_partition->head->next;
	while (node1 != information->partition_B->cells_in_partition->tail)
	{
		struct cell *cell1 = (struct cell *)node1->data_structure;
		area1 = 0;
		if (cell1->which_partition == 1)
			area1 = cell1->area2;
		printf("%d      %d           %4d    %4d\n", cell1->identifier + 1, cell1->which_partition, cell1->area, cell1->area2); //
		initial_partition_areas_B = initial_partition_areas_B + area1;
		node1 = node1->next;
	}
	printf("partition 0 area: %4d \npartition 1 area: %4d \n", initial_partition_areas_A, initial_partition_areas_B); //////////////////初始分區面積

	/////////////////
	///////////////
	// FILE *data = fopen(information->results_filename, "a");
	while (cells_can_still_be_moved)
	{

		if (PRINT_PARTITION_STATES || RUN_DEMO_WITH_TESTDATA)
		{
			// print_gain_arrays(information->access_[PARTITION_A]);//暫時關閉
			// print_gain_arrays(information->access_[PARTITION_B]);//暫時關閉
			//////////////
			/*print_partition_cell_ids(information->partition_A);
			print_partition_cell_ids(information->partition_B);*/
			//////////////
		}
		// fprintf(data, "%d, %d\n", timestep, information->current_cutstate);
		cells_can_still_be_moved = FM_pass(information);
		/////////////////////////////
		/*
		printf("partition %d 's Cell ID: ", information->partition_B->which_partition);
		struct node *current_node = information->partition_B->cells_in_partition->head->next;
		while (current_node != information->partition_B->cells_in_partition->tail)
		{
			struct cell *cell = (struct cell *)current_node->data_structure;
			printf("%d ", cell->identifier); //
			current_node = current_node->next;
		}
		printf("\n");

		/////////////////////////////
		printf("partition %d 's Cell ID: ", information->partition_A->which_partition);
		struct node *current_node1 = information->partition_A->cells_in_partition->head->next;
		while (current_node1 != information->partition_A->cells_in_partition->tail)
		{
			struct cell *cell = (struct cell *)current_node1->data_structure;
			printf("%d ", cell->identifier); //
			current_node1 = current_node1->next;
		}
		printf("\n");*/
		/*print_partition_cell_ids(information->partition_A);
		print_partition_cell_ids(information->partition_B);*/
		/////////////////////////////
		if (DOUBLE_CHECK_CUTSTATE_VALUES)
			check_cutstate_values(information);

		if (FM_CUTOFF && (information->current_cutstate > (information->lowest_cutstate + FM_CUTOFF_THRESHOLD)))
			break;

		timestep++;
	}

	// fclose(data);
	if (FM_REPEAT)
	{
		int i;
		for (i = 0; i < information->CELL_array_size; i++)
		{
			information->FM_chromosome->gene_array[i] = information->CELL_array[i]->which_partition;
		}
	}
}

// Go through each net in NET_array, check to see if net has at least one cell in each
// Because this goes through the entire netlist at every pass, it should not be used during normal cycles
void check_cutstate_values(struct condensed *information)
{

	struct net *temp_net;
	int cutstate_count = 0;
	int i;
	for (i = 0; i < information->NET_array_size; i++)
	{
		temp_net = information->NET_array[i];
		if ((temp_net->num_cells_in_[PARTITION_A] > 0) && (temp_net->num_cells_in_[PARTITION_B] > 0))
			cutstate_count++;
	}
	printf("Checked cutstate value: %d\n", cutstate_count);
}

// Prints the identifers of the cells in each list
void print_gain_arrays(struct partition *partition)
{
	int i;
	printf("*********\n");
	for (i = 0; i < partition->GAIN_array_size; i++)
	{
		printf("- ");
		print_dll(partition->GAIN_array[i], CELL);
	}
}

// Calls calculate_initial_gains for both partitions
void calculate_initial_gains_wrapper(struct condensed *information)
{
	calculate_initial_gains(information->partition_A, PARTITION_A, information->max_nets, information->FM_chromosome);
	calculate_initial_gains(information->partition_B, PARTITION_B, information->max_nets, information->FM_chromosome);
}

void calculate_initial_gains(struct partition *partition, partition_type label, int max_nets, struct chromosome *FM_chromosome)
{
	// Access variables
	struct node *temp_cell_node;
	struct node *temp_net_node;
	struct dll *temp_netlist;
	struct cell *temp_cell;
	struct net *temp_net;
	// Placement variables
	int cell_placement;
	struct dll **GAIN_array = partition->GAIN_array;
	struct dll *cell_gain_dll;

	// Access the first node in the partition with cell data
	temp_cell_node = partition->cells_in_partition->head->next;
	// Go through each cell, calculate gain, add to correct place in GAIN_array
	while (temp_cell_node != partition->cells_in_partition->tail)
	{
		temp_cell = (struct cell *)(temp_cell_node->data_structure);
		temp_netlist = temp_cell->nets;
		temp_net_node = ((struct node *)(temp_netlist->head))->next;
		while (temp_net_node != temp_netlist->tail)
		{
			temp_net = temp_net_node->data_structure;
			// Check if the cell is on a critical net, change gain on whether it would be good/bad to move cell
			if (temp_net->num_cells_in_[label] == 1)
			{
				temp_cell->gain += 1;
			}
			if (temp_net->num_cells_in_[!label] == 0)
			{
				temp_cell->gain -= 1;
			}
			temp_net_node = temp_net_node->next;
		}

		// Now that initial gain has been calculated, add to the appropriate dll in GAIN_array
		cell_placement = max_nets + temp_cell->gain;
		cell_gain_dll = GAIN_array[cell_placement];
		temp_cell->GAIN_array_node = insert_node(cell_gain_dll, 0, temp_cell);
		temp_cell_node = temp_cell_node->next;
	}

	// Set max_gain for partition
	update_max_gain_pointer(partition);
}

void update_max_gain_pointer(struct partition *partition)
{
	// Set max_gain
	int i;
	int no_pointer_change = 1;
	struct dll *temp_gain_list;
	// Go through the highest gains first, check dll sizes until the first nonzero gain list
	for (i = (partition->GAIN_array_size - 1); i >= 0; i--)
	{
		temp_gain_list = partition->GAIN_array[i];
		if (temp_gain_list->size < 0)
		{
			printf("Gain array list LESS THAN zero!\n");
		}
		if (temp_gain_list->size > 0)
		{
			partition->max_gain_pointer = ((struct node *)(temp_gain_list->head))->next;
			no_pointer_change = 0;
			break;
		}
	}
	// If no more cells in the partition, set the pointer to NULL
	if (no_pointer_change)
	{
		partition->max_gain_pointer = NULL;
	}
}

// #####################################################################################################

// When updating the gains, remember to check for nets that are added/removed from cutstate
// Each pass goes through one cell's netlist 4 times, checking for different conditions to update the cell's gains
// This function choses one cell to move to the other partition, making sure that balance is not upset.
//   It then updates the gains of cells on connected nets
int round = 1; /////////////////分割圈數
int FM_pass(struct condensed *information)
{

	// Access variables;
	struct dll *base_cell_netlist;
	struct node *temp_net_node;
	struct net *temp_net;
	struct node *temp_cell_node;

	// The cutsize will probably change after this code (hopefully to a lower number)
	int new_cutsize = information->current_cutstate;

	// base_cell is the notation for the cell to be moved between partitions
	struct cell *base_cell;
	struct node *node_to_be_freed;

	// Check both partitions with the max_gain_pointer with the higher cell
	struct node *node_A = information->partition_A->max_gain_pointer;
	struct node *node_B = information->partition_B->max_gain_pointer;

	// If no more cells, the algorithm halts
	if (node_A == NULL && node_B == NULL)
	{
		return 0;
	}

	struct cell *cell_A;
	struct cell *cell_B;

	// Make sure null pointers aren't dereferenced
	if (node_A != NULL)
	{
		cell_A = node_A->data_structure;
	}

	if (node_B != NULL)
	{
		cell_B = node_B->data_structure;
	}

	// Chose base_cell
	// Either other base_cell is not a candidate or (cell has the highest gain and does not mess up balance)
	// printf("partition A0 area: %4d \npartition B1 area: %4d \n", information->partition_A->total_partition_area, information->partition_B->total_partition_area); ////面積測試
	if ((node_A != NULL) && (information->partition_A->total_partition_area - cell_A->area > /*information->desired_area*/ ((information->partition_A->total_partition_area + information->partition_B->total_partition_area) / 2) - information->tolerance))
	{
		base_cell = cell_A;
		node_to_be_freed = node_A;
	}
	else if ((node_B != NULL) && (information->partition_B->total_partition_area - cell_B->area2 /*改成第二製程面積*/ > /*information->total_area - information->desired_area*/ ((information->partition_A->total_partition_area + information->partition_B->total_partition_area) / 2) - information->tolerance1))
	{
		base_cell = cell_B;
		node_to_be_freed = node_B;
	}
	else
	{
		printf("No cells can be chosen\n");
		return 0;
	}

	// Set origin and destination references
	int base_cell_origin, base_cell_destination;

	base_cell_origin = base_cell->which_partition;
	base_cell_destination = !(base_cell_origin);

	// Switch the base_cell partition
	base_cell->which_partition = base_cell_destination;
	base_cell->partition = information->access_[base_cell_destination];
	printf("chosen cell : %d \n", base_cell->identifier + 1);  /////////////
	printf("switch to partition %d\n", base_cell_destination); ///////////
	/*	FILE *shmetisInput = fopen("Output_minimumcut", "a");	   //////////////寫入Output
	fprintf(shmetisInput, "cell   partition\n");			   /////////////寫入Output
	fprintf(shmetisInput, "cut: %d\n", new_cutsize);		   /////////////寫入Output
	printf("cell   partition\n");*/
	/////////////////////////////
	// printf("partition %d 's Cell ID: \n", information->partition_A->which_partition);
	/*	struct node *current_node = information->partition_A->cells_in_partition->head->next;
		while (current_node != information->partition_A->cells_in_partition->tail)
		{
			struct cell *cell = (struct cell *)current_node->data_structure;
			// if (cell->which_partition == 0)
			//{
			printf("%d      %d \n", cell->identifier, cell->which_partition);				 //
			fprintf(shmetisInput, "%d      %d \n", cell->identifier, cell->which_partition); //////////寫入Output
			//}
			current_node = current_node->next;
		}*/
	//////////////////
	/*if (base_cell_destination == 0)
	{
		printf("%d ", base_cell->identifier);
	}*/
	// printf("\n");
	/*print_partition_cell_ids(information->partition_A);
	print_partition_cell_ids(information->partition_B);*/
	/////////////////////////////
	// printf("partition %d 's Cell ID: \n", information->partition_B->which_partition);
	/*	struct node *current_node1 = information->partition_B->cells_in_partition->head->next;
		while (current_node1 != information->partition_B->cells_in_partition->tail)
		{
			struct cell *cell = (struct cell *)current_node1->data_structure;
			// if (cell->which_partition == 1)
			//{
			printf("%d      %d \n", cell->identifier, cell->which_partition);				 //
			fprintf(shmetisInput, "%d      %d \n", cell->identifier, cell->which_partition); //////////寫入Output
			//}
			current_node1 = current_node1->next;
		}*/
	/////////////
	// fclose(shmetisInput); // 寫入Output
	/////////////
	/*if (base_cell_destination == 1)
	{
		printf("%d ", base_cell->identifier);
	}*/
	printf("\n");
	/*print_partition_cell_ids(information->partition_A);
	print_partition_cell_ids(information->partition_B);*/
	/////////////////////////////
	// update FM_chromosome (for Genetic Algorithm)
	int *gene_array;
	if (FM_REPEAT)
	{
		gene_array = information->FM_chromosome->gene_array;
		gene_array[base_cell->identifier] = base_cell_destination;
	}

	// update partitions now that base_cell is chosen
	// Find the dll the base_cell belongs to (GAIN_array_size is the same for both partitions)
	int position = information->access_[base_cell_origin]->GAIN_array_size / 2 + base_cell->gain;
	// Update partition areas
	information->access_[base_cell_origin]->total_partition_area -= base_cell->area;
	information->access_[base_cell_destination]->total_partition_area += base_cell->area;
	// Delete node from Gain array
	free(remove_node(node_to_be_freed, information->access_[base_cell_origin]->GAIN_array[position]));
	base_cell->GAIN_array_node = NULL;
	base_cell->cell_state = LOCKED;
	// Setup for the while loop
	base_cell_netlist = base_cell->nets;
	temp_net_node = base_cell_netlist->head->next;

	int num_cells_in_destination, num_cells_in_origin;

	// For each net in the base cell netlist
	while (temp_net_node != base_cell_netlist->tail)
	{
		// Access the net
		temp_net = temp_net_node->data_structure;

		// Skip nets with only one cell
		if (temp_net->free_cells->size == 1)
		{
			temp_net_node = temp_net_node->next;
			continue;
		}

		// if the net was previously not in the cutstate, increment the gains of all cells in net, increase cutsize
		if (temp_net->num_cells_in_[base_cell_destination] == 0)
		{
			new_cutsize = change_gain_of_cell_in_net(base_cell, ALL, INCREMENT, new_cutsize, temp_net->free_cells, base_cell_destination, temp_net->locked_cells, information);
		}
		// else if there was one cell previously in destination partition, decrement the gain of that cell
		else if (temp_net->num_cells_in_[base_cell_destination] == 1)
		{
			new_cutsize = change_gain_of_cell_in_net(base_cell, ONE, DECREMENT, new_cutsize, temp_net->free_cells, base_cell_destination, temp_net->locked_cells, information);
		}

		// Update net count now that base_cell has switched partitions
		temp_net->num_cells_in_[base_cell_destination] += 1;
		temp_net->num_cells_in_[base_cell_origin] -= 1;

		// if the origin partition now has no cells, decrement the gains of all cells in net, decrease cutsize
		if (temp_net->num_cells_in_[base_cell_origin] == 0)
		{
			new_cutsize = change_gain_of_cell_in_net(base_cell, ALL, DECREMENT, new_cutsize, temp_net->free_cells, base_cell_destination, temp_net->locked_cells, information);
		}
		else if (temp_net->num_cells_in_[base_cell_origin] == 1)
		{
			new_cutsize = change_gain_of_cell_in_net(base_cell, ONE, INCREMENT, new_cutsize, temp_net->free_cells, base_cell_destination, temp_net->locked_cells, information);
		}

		temp_net_node = temp_net_node->next;
	}

	// Print current cutstate
	if (PRINT_PASS_CUTSTATE_VALUES || RUN_DEMO_WITH_TESTDATA)
	{
		printf("Pass cutstate value: %d\n", new_cutsize);
	}

	// Check and update cutstate value
	if (new_cutsize < information->lowest_cutstate)
	{
		information->lowest_cutstate = new_cutsize;
	}
	information->current_cutstate = new_cutsize;
	////////////
	FILE *shmetisInput = fopen("Output_minimumcut", "a"); //////////////寫入Output
	if (new_cutsize <= information->lowest_cutstate)
	{
		fprintf(shmetisInput, "cell   partition\n");					  /////////////寫入Output
		fprintf(shmetisInput, "cut: %d  Round:%d\n", new_cutsize, round); /////////////寫入Output
	}
	round++; ////圈數加1
	printf("cell   partition   area    !area    gain\n");
	long total_area_in_partition_A = 0; /////////////
	long total_area_in_partition_B = 0; /////////////
	long total_area_in_partition_C = 0; /////////////
	long total_area_in_partition_D = 0; /////////////
	int area;
	int area2;
	int area3;
	int area4;
	// printf("partition %d 's Cell ID: \n", information->partition_A->which_partition);
	struct node *current_node = information->partition_A->cells_in_partition->head->next;
	while (current_node != information->partition_A->cells_in_partition->tail)
	{
		///////
		area = 0;
		area2 = 0;
		area3 = 0;
		area4 = 0;
		///////
		struct cell *cell = (struct cell *)current_node->data_structure;
		if (cell->which_partition == 0)
		{
			area = cell->area;
			area3 = cell->area2;
		}
		else
		{
			area2 = cell->area2;
			area4 = cell->area;
		}
		printf("%d      %d           %4d    %4d    %2d\n", cell->identifier + 1, cell->which_partition, cell->area, cell->area2, cell->gain); //
		total_area_in_partition_A = total_area_in_partition_A + area;
		total_area_in_partition_C = total_area_in_partition_C + area3;
		total_area_in_partition_B = total_area_in_partition_B + area2;
		total_area_in_partition_D = total_area_in_partition_D + area4;
		if (new_cutsize <= information->lowest_cutstate)
		{
			fprintf(shmetisInput, "a%d      %d \n", cell->identifier + 1, cell->which_partition); //////////寫入Output
		}
		//}
		current_node = current_node->next;
	}
	struct node *current_node1 = information->partition_B->cells_in_partition->head->next;
	while (current_node1 != information->partition_B->cells_in_partition->tail)
	{
		///////
		area = 0;
		area2 = 0;
		area3 = 0;
		area4 = 0;
		///////
		struct cell *cell2 = (struct cell *)current_node1->data_structure;
		if (cell2->which_partition == 1)
		{
			area2 = cell2->area2;
			area4 = cell2->area;
		}
		else
		{
			area = cell2->area;
			area3 = cell2->area2;
		}
		printf("%d      %d           %4d    %4d    %2d\n", cell2->identifier + 1, cell2->which_partition, cell2->area, cell2->area2, cell2->gain); //
		total_area_in_partition_A = total_area_in_partition_A + area;
		total_area_in_partition_C = total_area_in_partition_C + area3;
		total_area_in_partition_B = total_area_in_partition_B + area2;
		total_area_in_partition_D = total_area_in_partition_D + area4;
		if (new_cutsize <= information->lowest_cutstate)
		{
			fprintf(shmetisInput, "a%d      %d \n", cell2->identifier + 1, cell2->which_partition); //////////寫入Output
		}
		//}
		current_node1 = current_node1->next;
	}
	if (new_cutsize <= information->lowest_cutstate)
	{
		fprintf(shmetisInput, "A: %4ld or %4ld\nB: %4ld or %4ld\n", total_area_in_partition_A, total_area_in_partition_C, total_area_in_partition_B, total_area_in_partition_D); //////////寫入Output
	}
	fclose(shmetisInput); // 寫入Output
	printf("A: %4ld or %4ld\n", total_area_in_partition_A, total_area_in_partition_C);
	printf("B: %4ld or %4ld\n", total_area_in_partition_B, total_area_in_partition_D);

	////////////
	information->partition_A->total_partition_area = total_area_in_partition_A; // 面積重新賦值
	information->partition_B->total_partition_area = total_area_in_partition_B; // 面積重新賦值
	//  Set the max_gain pointers to the new highest cells
	update_max_gain_pointer(information->partition_A);
	update_max_gain_pointer(information->partition_B);

	return 1;
}

int change_gain_of_cell_in_net(struct cell *base_cell, cells_to_change scope, change_direction operation, int old_cutsize, struct dll *cellist, int partition_with_isolated_cell, struct dll *locked_cells, struct condensed *information)
{
	struct node *temp_cell_node = cellist->head->next;

	struct cell *temp_cell;
	int position, old_gain;
	partition_type temp_cell_origin;

	struct node *locked_cell_node;

	// Go through every cell in the net's cellist
	while (temp_cell_node != cellist->tail)
	{
		// Access the cell
		temp_cell = temp_cell_node->data_structure;
		temp_cell_origin = temp_cell->which_partition;

		// Remove locked cells from net's free_cell list
		if (temp_cell == base_cell || temp_cell->cell_state == LOCKED)
		{
			// Move to next node, get rid of bad one
			temp_cell_node = temp_cell_node->next;
			locked_cell_node = remove_node(temp_cell_node->previous, cellist);
			insert_node(locked_cells, 0, locked_cell_node->data_structure);
			free(locked_cell_node);
			continue;
		}

		old_gain = temp_cell->gain;

		if (scope == ONE)
		{
			// Update gain
			if (temp_cell->which_partition == partition_with_isolated_cell)
			{
				if (operation == INCREMENT)
					temp_cell->gain++;
				else
					temp_cell->gain--;
			}
		}
		else if (scope == ALL)
		{
			// Update gain
			if (operation == INCREMENT)
				temp_cell->gain++;
			else
				temp_cell->gain--;
		}

		position = information->access_[temp_cell_origin]->GAIN_array_size / 2 + old_gain;
		// Remove temp_cell from old position in GAIN_array
		free(remove_node(temp_cell->GAIN_array_node, information->access_[temp_cell_origin]->GAIN_array[position]));
		// Find the new position temp_cell should be in based on it's gain
		position = information->access_[temp_cell_origin]->GAIN_array_size / 2 + temp_cell->gain;
		// Insert the cell and store the node information
		temp_cell->GAIN_array_node = insert_node(information->access_[temp_cell_origin]->GAIN_array[position], 0, temp_cell);

		// Move to next cell
		temp_cell_node = temp_cell_node->next;
	}

	int new_cutsize = old_cutsize;
	// Update cutsize
	if (scope == ALL)
	{
		if (operation == INCREMENT)
			new_cutsize++;
		else
			new_cutsize--;
	}
	return new_cutsize;
}
