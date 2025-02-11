//a small helper function, does binary search for target on sorted array for part of the array given by start and end
int binary_search_part(int *arr, int start, int end, int target) {
    // TODO: replace with stdlib.h 
    //void* bsearch(const void* key, const void* base, size_t len,
    //size_t elem_size, <compare_function>);
    int low = start;
    int high = end - 1;

    while (low <= high) {
        int mid = low + (high - low) / 2;

        if (arr[mid] == target) {
            return mid; // Target found
        } else if (arr[mid] < target) {
            low = mid + 1; // Search upper half
        } else {
            high = mid - 1; // Search lower half
        }
    }

    return -1; // Target not found
}
