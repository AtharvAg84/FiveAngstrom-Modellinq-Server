import time
start_time = time.time()

from rcsbapi.search import TextQuery

# Search for structures associated with the phrase "Hemoglobin"

query = TextQuery(value="Aspirin    ")

# Execute the query by running it as a function
results = query()

# Results are returned as an iterator of result identifiers.
for rid in results:
    print(rid)
end_time = time.time()
print(f"Program started at: {time.ctime(start_time)}")
print(f"Program ended at: {time.ctime(end_time)}")
print(f"Total execution time: {end_time - start_time:.2f} seconds")

