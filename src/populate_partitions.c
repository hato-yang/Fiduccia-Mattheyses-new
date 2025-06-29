#include "../include/main.h"

// Assumes no particular order for the CELL_array
// This seems to work well when the majority of cells are roughly the same, with large outliers
void segregate_cells_randomly(struct condensed *information)
{
	// Cells are divided randomly into two temporary partiton lists
	//	ex. If rand() >0.5, cell(i) goes to A, otherwise goes to B
	//  The cell areas are meanwhile added to the two totals
	//  At the end, the total areas are compared. If they don't match the desired ratio,
	//	lists are deleted and computer tries again.
	// If desired ratio met, the temporary lists are transcribed into the partition data_structures

	int i;

	struct dll *list_of_cells_A;
	struct dll *list_of_cells_B;
	/////////////////計算分布機率
	float x;
	x = (float)information->total_area / (information->total_area + information->total_area2);
	/////////////////
	/*int*/ long total_partition_areas[2];
	total_partition_areas[PARTITION_A] = 0;
	total_partition_areas[PARTITION_B] = 0;

	struct cell **CELL_array = information->CELL_array;

	// which partition the cell is assigned to PARTITION_A, PARTITION_B
	int partition_placement;
	// Loop until acceptable balance is found
	while (1)
	{
		list_of_cells_A = malloc(sizeof(struct dll));
		list_of_cells_B = malloc(sizeof(struct dll));

		// Seed random for different results (adding rand() prevents the seed from being reused in short time intervals)
		// Initialize the dlls
		initialize_dll(list_of_cells_A);
		initialize_dll(list_of_cells_B);

		srand(time(NULL) + rand());
		//////////////////處理順序洗牌
		// 建立亂數 index array
		int *index_array = malloc(sizeof(int) * information->CELL_array_size);
		for (int i = 0; i < information->CELL_array_size; i++)
		{
			index_array[i] = i;
		}
		// Fisher–Yates 洗牌
		for (int i = information->CELL_array_size - 1; i > 0; i--)
		{
			int j = rand() % (i + 1);
			int temp = index_array[i];
			index_array[i] = index_array[j];
			index_array[j] = temp;
		}
		//////////////////
		// Assign every cell to a partition
		for (i = 0; i < information->CELL_array_size; i++)
		{
			//////////// 用隨機的 index 處理 cell
			int idx = index_array[i];
			struct cell *cell = CELL_array[idx];
			////////////
			// Decide which partition the cell will go into
			// First option is random, repeat runs are based on FM_chromosome
			if (information->FM_chromosome != NULL)
			{
				partition_placement = information->FM_chromosome->gene_array[i];
			}
			else
			{
				// float random_value = (float)rand() / RAND_MAX; // 生成 0 到 1 的隨機數
				// partition_placement = (random_value < x);	   // 如果隨機數小於 x，則選擇 partition A
				//     根據當前分區面積動態決定分配到哪個分區
				if (/*total_partition_areas[PARTITION_A] + CELL_array[i] /*改成 cell*->area <= total_partition_areas[PARTITION_B] + CELL_array[i] /*改成 cell*->area2 &&*/ total_partition_areas[PARTITION_A] + CELL_array[i] /*改成 cell*/->area < ((float)information->utilA / 100 * information->die_area))
				{
					/*if (total_partition_areas[PARTITION_A] + CELL_array[i] /*改成 cell->area >= ((float)information->utilA / 200 * information->die_area) && total_partition_areas[PARTITION_B] + CELL_array[i] /*改成 cell->area2 < ((float)information->utilB / 100 * information->die_area))
					{
						partition_placement = PARTITION_B;
					}
					else
					{*/
					partition_placement = PARTITION_A;
					//}
				}
				else
				{
					partition_placement = PARTITION_B;
				}
			}
			// partition_placement = ((rand() % (int)(1.0 / (information->ratio))) == 0); // 原始隨機分配cell方法

			// total_partition_areas[partition_placement] += CELL_array[i]->area; // 原始代入面積寫法
			////////////////
			if (partition_placement == PARTITION_A)
			{
				total_partition_areas[partition_placement] += CELL_array[i] /*改成 cell*/->area;
				// printf("cell# : %d cellarea : %d total.A : %ld\n", i, CELL_array[i]->area, total_partition_areas[partition_placement]);
			}
			if (partition_placement == PARTITION_B)
			{
				total_partition_areas[partition_placement] += CELL_array[i] /*改成 cell*/->area2;
				// printf("cell# : %d cellarea : %d total.B : %ld\n", i, CELL_array[i]->area2, total_partition_areas[partition_placement]);
			}
			////////////////
			CELL_array[i] /*改成 cell*/->partition = information->access_[partition_placement];
			CELL_array[i] /*改成 cell*/->which_partition = partition_placement;

			if (partition_placement == PARTITION_A)
			{
				insert_node(list_of_cells_A, 0, CELL_array[i] /*改成 cell*/);
			}
			else
			{
				insert_node(list_of_cells_B, 0, CELL_array[i] /*改成 cell*/);
			}
			// printf("CELL_array[%d]: %d partition: %d while partition: %d\n", i, CELL_array[i]->identifier + 1, CELL_array[i]->partition->which_partition, CELL_array[i]->which_partition);
		}
		/*int*/ long two_partition_area = 0;														  // 初始總面積計算
		two_partition_area = total_partition_areas[PARTITION_A] + total_partition_areas[PARTITION_B]; // 初始總面積計算
		// printf("%ld   %ld   %ld   %ld\n", total_partition_areas[PARTITION_A], total_partition_areas[PARTITION_B], two_partition_area, two_partition_area / 2);
		// If the partition is within tolerance, break the loop and save to partition structs
		// Otherwise free dlls and try again.
		// if (total_partition_areas[PARTITION_A] < ((two_partition_area / 2) + (information->tolerance_Macro * 3)) && total_partition_areas[PARTITION_A] > ((two_partition_area / 2) - (information->tolerance_Macro * 3)))
		if (total_partition_areas[PARTITION_A] < ((double)information->utilA / 100 * information->die_area + (information->tolerance_Macro * 2)) && total_partition_areas[PARTITION_B] < ((double)information->utilB / 100 * information->die_area + (information->tolerance1_Macro * 2)))
		// if (total_partition_areas[PARTITION_A] < 2838171970 && total_partition_areas[PARTITION_B] < 2838171970)//兩邊面積皆小於Die size
		{
			// printf("%ld,%f,%d,%f,%d\n", information->die_area, (float)information->utilA / 100 * information->die_area, (information->tolerance_Macro), (float)information->utilB / 100 * information->die_area, (information->tolerance1_Macro));
			information->initial_area = total_partition_areas[PARTITION_A];
			information->initial_area2 = total_partition_areas[PARTITION_B];
			break;
		}
		// 原始代入面積寫法  information->desired_area 修改成 two_partition_area/2
		// Get ready for next loop
		garbage_collection_dll(list_of_cells_A, DO_NOT_DEALLOC_DATA);
		garbage_collection_dll(list_of_cells_B, DO_NOT_DEALLOC_DATA);
		total_partition_areas[PARTITION_A] = 0;
		total_partition_areas[PARTITION_B] = 0;
		free(index_array);
	}

	copy_cells_into_partitions(information->partition_A, information->partition_B, list_of_cells_A, list_of_cells_B, total_partition_areas[PARTITION_A], total_partition_areas[PARTITION_B]);

	garbage_collection_dll(list_of_cells_A, DO_NOT_DEALLOC_DATA);
	garbage_collection_dll(list_of_cells_B, DO_NOT_DEALLOC_DATA);
}

// Add cells, return dll of cutstate nets
void copy_cells_into_partitions(struct partition *partition_A, struct partition *partition_B, struct dll *list_of_cells_A, struct dll *list_of_cells_B, long total_partition_area_A, long total_partition_area_B)
{

	partition_A->total_partition_area = total_partition_area_A;

	// Go through each cell. Each cell has a list of nets. Go through each net. Add 1 to the number_cells_in_partition_X;
	struct node *temp_node;
	struct cell *temp_cell;

	// Access the first data node
	temp_node = ((struct node *)list_of_cells_A->head)->next;
	// Go through list of cells
	while (temp_node != list_of_cells_A->tail)
	{
		temp_cell = temp_node->data_structure;
		// Add cell to partition
		insert_node(partition_A->cells_in_partition, 0, temp_cell);
		update_net_partition_count(temp_cell, PARTITION_A);
		// Move to next cell
		temp_node = temp_node->next;
	}

	partition_B->total_partition_area = total_partition_area_B;

	// Access the first data node
	temp_node = ((struct node *)list_of_cells_B->head)->next;
	// Go through list of cells
	while (temp_node != list_of_cells_B->tail)
	{
		temp_cell = temp_node->data_structure;
		// Add cell to partition
		insert_node(partition_B->cells_in_partition, 0, temp_cell);
		update_net_partition_count(temp_cell, PARTITION_B);
		// Move to next cell
		temp_node = temp_node->next;
	}
}

// Go through each net and increment the appropriate counter variable
void update_net_partition_count(struct cell *assigned_cell, partition_type partition)
{
	struct dll *netlist = assigned_cell->nets;
	struct node *temp_net_node = ((struct node *)netlist->head)->next;
	struct net *temp_net;
	while (temp_net_node != netlist->tail)
	{
		temp_net = temp_net_node->data_structure;

		temp_net->num_cells_in_[partition] += 1;
		// Move to next net
		temp_net_node = temp_net_node->next;
	}
}

int calculate_initial_cutstate(struct net **NET_array, int NET_array_size, struct condensed *information)
{
	// Go through each net in NET_array, check to see if net has at least one cell in each

	struct net *temp_net;
	int cutstate_count = 0;
	int i;
	for (i = 0; i < NET_array_size; i++)
	{
		temp_net = NET_array[i];
		// printf("Net %d: A=%d B=%d\n", i, temp_net->num_cells_in_[0], temp_net->num_cells_in_[1]); /// 計算net在兩partition中個別的cell數
		if ((temp_net->num_cells_in_[PARTITION_A] > 0) && (temp_net->num_cells_in_[PARTITION_B] > 0))
			cutstate_count++;
	}
	printf("Initial cutstate value: %d\n", cutstate_count);
	information->current_cutstate = cutstate_count;
	return cutstate_count;
}

void segregate_cells_by_net_order(struct condensed *info)
{

	for (int i = 0; i < info->CELL_array_size; i++)
	{
		info->CELL_array[i]->which_partition = -1;
		info->CELL_array[i]->partition = NULL;
	}
	int i;

	struct dll *list_of_cells_A;
	struct dll *list_of_cells_B;

	long total_partition_areas[2];
	total_partition_areas[PARTITION_A] = 0;
	total_partition_areas[PARTITION_B] = 0;
	int partition_placement;

	while (1)
	{
		list_of_cells_A = malloc(sizeof(struct dll));
		list_of_cells_B = malloc(sizeof(struct dll));

		initialize_dll(list_of_cells_A);
		initialize_dll(list_of_cells_B);

		int threshold = info->max_cell_count * 0.8;
		int large_net_count = 0;
		for (int i = 0; i < info->NET_array_size; i++)
		{
			struct net *n = info->NET_array[i];
			struct node *cur = n->free_cells->head->next;

			if (n->number_of_cells <= threshold)
			{
				int target_partition = rand() % 2;
				while (cur != n->free_cells->tail)
				{
					struct cell *c = (struct cell *)cur->data_structure;
					// printf("SMALL NET: Net %d, Cell %d, Assigned? %d, Area: %d, Area2: %d\n", n->identifier, c->identifier + 1, c->which_partition, c->area, c->area2);

					if (c->which_partition == -1)
					{
						float maxA = (float)info->utilA / 100 * info->die_area;
						float maxB = (float)info->utilB / 100 * info->die_area;

						if (target_partition == 0)
						{
							if (total_partition_areas[PARTITION_A] + c->area <= maxA)
							{
								c->which_partition = target_partition;
								c->partition = info->access_[target_partition];
								insert_node(list_of_cells_A, 0, c);
								total_partition_areas[PARTITION_A] += c->area;
								// printf("  → Assigned to %d, new total: %ld\n", target_partition, total_partition_areas[PARTITION_A]);
							}
							else
							{
								c->which_partition = !target_partition;
								c->partition = info->access_[!target_partition];
								insert_node(list_of_cells_B, 0, c);
								total_partition_areas[PARTITION_B] += c->area2;
								// printf("  → Assigned to %d, new total: %ld\n", !target_partition, total_partition_areas[PARTITION_B]);
							}
						}
						else
						{
							c->which_partition = target_partition;
							c->partition = info->access_[target_partition];
							insert_node(list_of_cells_B, 0, c);
							total_partition_areas[PARTITION_B] += c->area2;
							// printf("  → Assigned to %d, new total: %ld\n", target_partition, total_partition_areas[PARTITION_B]);
						}
					}
					cur = cur->next;
				}
			}
			else
			{
				large_net_count++;
				int a_count = 0, b_count = 0;
				cur = n->free_cells->head->next;
				while (cur != n->free_cells->tail)
				{
					struct cell *c = (struct cell *)cur->data_structure;
					// printf("LARGE NET: Net %d, Cell %d, Assigned? %d\n", n->identifier, c->identifier + 1, c->which_partition);

					if (c->which_partition == -1)
					{
						int assign_to = (a_count <= b_count) ? 0 : 1;
						c->which_partition = assign_to;
						c->partition = info->access_[assign_to];

						if (assign_to == 0)
						{
							a_count++;
							insert_node(list_of_cells_A, 0, c);
							total_partition_areas[PARTITION_A] += c->area;
							// printf("  → Assigned to %d (large net), new total: %ld\n", assign_to, total_partition_areas[PARTITION_A]);
						}
						else
						{
							b_count++;
							insert_node(list_of_cells_B, 0, c);
							total_partition_areas[PARTITION_B] += c->area2;
							// printf("  → Assigned to %d (large net), new total: %ld\n", assign_to, total_partition_areas[PARTITION_B]);
						}
					}
					cur = cur->next;
				}
			}
		}

		double targetA = (double)info->die_area * info->utilA / 100.0 + info->tolerance_Macro * 2;
		double targetB = (double)info->die_area * info->utilB / 100.0 + info->tolerance1_Macro * 2;

		// printf("CHECK: areaA = %ld (limit %.2f), areaB = %ld (limit %.2f)\n",total_partition_areas[PARTITION_A], targetA,total_partition_areas[PARTITION_B], targetB);

		if (total_partition_areas[PARTITION_A] < targetA && total_partition_areas[PARTITION_B] < targetB)
		{
			printf("#max_nets_cell: %d, threshold: %d, #large net: %d\n", info->max_cell_count, threshold, large_net_count);
			printf("==> Valid partition found.\n");
			info->initial_area = total_partition_areas[PARTITION_A];
			info->initial_area2 = total_partition_areas[PARTITION_B];
			break;
		}

		garbage_collection_dll(list_of_cells_A, DO_NOT_DEALLOC_DATA);
		garbage_collection_dll(list_of_cells_B, DO_NOT_DEALLOC_DATA);
		total_partition_areas[PARTITION_A] = 0;
		total_partition_areas[PARTITION_B] = 0;
	}

	copy_cells_into_partitions(info->partition_A, info->partition_B, list_of_cells_A, list_of_cells_B, total_partition_areas[PARTITION_A], total_partition_areas[PARTITION_B]);

	garbage_collection_dll(list_of_cells_A, DO_NOT_DEALLOC_DATA);
	garbage_collection_dll(list_of_cells_B, DO_NOT_DEALLOC_DATA);
}
