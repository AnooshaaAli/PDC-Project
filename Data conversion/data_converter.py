from collections import defaultdict
import random

# Step 1: Load the social network edges
edge_activities = defaultdict(lambda: {'retweet': 0, 'reply': 0, 'mention': 0})

with open('social_network.edgelist', 'r') as f:
    for line in f:
        u, v = map(int, line.strip().split())
        edge_activities[(u, v)]  # Initialize with zeros

# Step 2: Load the activity log and update counts
activity_map = {
    'RT': 'retweet',
    'RE': 'reply',
    'MT': 'mention'
}

with open('activity_time.txt', 'r') as f:
    for line in f:
        from_id, to_id, _, activity = line.strip().split()
        from_id = int(from_id)
        to_id = int(to_id)
        if (from_id, to_id) in edge_activities:
            edge_activities[(from_id, to_id)][activity_map[activity]] += 1

# Step 3: Write to output.txt (no header)
with open('output.txt', 'w') as f:
    for (u, v), counts in edge_activities.items():
        f.write(f"{u} {v} {counts['retweet']} {counts['reply']} {counts['mention']}\n")



# Load Gnutella edges
edges = []
with open('Gnutella.txt', 'r') as f:
    for line in f:
        u, v = map(int, line.strip().split())
        edges.append((u, v))

# Randomly assign activities
def random_activity():
    # 50% chance of all zeros, 30% chance of one activity, 20% chance of multiple
    prob = random.random()
    if prob < 0.5:
        return 0, 0, 0
    elif prob < 0.8:
        choice = random.choice(['retweet', 'reply', 'mention'])
        return (
            random.randint(1, 3) if choice == 'retweet' else 0,
            random.randint(1, 3) if choice == 'reply' else 0,
            random.randint(1, 3) if choice == 'mention' else 0
        )
    else:
        return random.randint(0, 3), random.randint(0, 3), random.randint(0, 3)

# Write to file
with open('dataset2.txt', 'w') as f:
    for u, v in edges:
        rt, re, mt = random_activity()
        f.write(f"{u} {v} {rt} {re} {mt}\n")

