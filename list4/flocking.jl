using Agents, LinearAlgebra
using DrWatson: @dict

@agent Bird ContinuousAgent{2} begin
    speed::Float64
    cohere_factor::Float64
    separation::Float64
    separate_factor::Float64
    match_factor::Float64
    visual_distance::Float64
    vision_angle::Float64
end

function initialize_model(;
    n_birds = 100,
    speed = 1.0,
    cohere_factor = 0.25,
    separation = 4.0,
    separate_factor = 0.25,
    match_factor = 0.01,
    visual_distance = 5.0,
    vision_angle = π/8,
    extent = (100, 100),
)
    space2d = ContinuousSpace(extent; spacing = visual_distance/1.5)
    model = ABM(Bird, space2d, scheduler = Schedulers.Randomly())
    for _ in 1:n_birds
        vel = Tuple(rand(model.rng, 2) * 2 .- 1)
        add_agent!(
            model,
            vel,
            speed,
            cohere_factor,
            separation,
            separate_factor,
            match_factor,
            visual_distance,
            vision_angle,
        )
    end
    return model
end

function angle_between(bird1, bird2)
    # This is the angle between the bird's heading and the vector between the birds
    heading1 = bird1.vel
    heading2 = bird2.pos .- bird1.pos
    heading1 = heading1 ./ norm(heading1)
    heading2 = heading2 ./ norm(heading2)
    dot_prod = dot(heading1, heading2)
    dot_prod = max(-1.0, min(1.0, dot_prod))
    angle = acos(dot_prod)
    return angle
end


function agent_step!(bird, model)
    # Obtain the ids of neighbors within the bird's visual distance
    all_neighbor_ids = nearby_ids(bird, model, bird.visual_distance)

    # neighbor_ids = all_neighbor_ids
    neighbor_ids = Int[]
    
    # clear the list of neighbors
    
    for id in all_neighbor_ids
        neighbor = model[id]
        # Only consider neighbors within the bird's visual angle
        angle = angle_between(bird, neighbor)
        if angle < bird.vision_angle
            push!(neighbor_ids, id)
        end
    end

    N = 0
    match = separate = cohere = (0.0, 0.0)
    # Calculate behaviour properties based on neighbors
    for id in neighbor_ids
        N += 1
        neighbor = model[id].pos
        heading = neighbor .- bird.pos

        # `cohere` computes the average position of neighboring birds
        cohere = cohere .+ heading
        if euclidean_distance(bird.pos, neighbor, model) < bird.separation
            # `separate` repels the bird away from neighboring birds
            separate = separate .- heading
        end
        # `match` computes the average trajectory of neighboring birds
        match = match .+ model[id].vel
    end
    N = max(N, 1)
    # Normalise results based on model input and neighbor count
    cohere = cohere ./ N .* bird.cohere_factor
    separate = separate ./ N .* bird.separate_factor
    match = match ./ N .* bird.match_factor
    # Compute velocity based on rules defined above
    bird.vel = (bird.vel .+ cohere .+ separate .+ match) ./ 2
    bird.vel = bird.vel ./ norm(bird.vel)
    # Move bird according to new velocity and speed
    move_agent!(bird, model, bird.speed)
end

# make the animation

using InteractiveDynamics
using GLMakie


function bird_marker(b::Bird)
    # x = b.visual_distance * cos(b.vision_angle)
    # y = b.visual_distance * sin(b.vision_angle)
    bird_polygon = Polygon(Point2f[(-0.5, -0.5), (1, 0), (-0.5, 0.5)])
    φ = atan(b.vel[2], b.vel[1]) #+ π/2 + π
    scale(rotate2D(bird_polygon, φ), 2)
    
end

model = initialize_model(
    n_birds = 300,
    speed = 1.0,
    extent = (200, 200),
    separation = 5.0,
    visual_distance = 5.0,
    vision_angle = π/4,
    cohere_factor = 0.25,
    separate_factor = 0.25,
    match_factor = 0.01,
)


figure, = abmplot(model; am = bird_marker, add_controls = true, title = "Flocking",
agent_step! = agent_step!)
figure

# abmvideo(
#     "flocking.mp4", model, agent_step!;
#     am = bird_marker,
#     framerate = 20, frames = 100,
#     title = "Flocking"
# )